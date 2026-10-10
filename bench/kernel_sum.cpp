// Micro-benchmark of the OPES kernel sum (getProbAndDerivatives without NLIST, PR cutoff scheme).
// Real kernels (a STATE file) are evaluated along a real CV path (a COLVAR file); all CVs are
// taken as periodic in [-pi,pi), as the torsions of ala2 and ala4.
//
// usage: kernel_sum STATE COLVAR [repeats]
// Variants:
//   current      the v2.9 code: vector of kernels, each with its own center and sigma vectors,
//                Value::difference, division by sigma, early exit after each CV
//   inv_sigma    same, multiplying by a stored 1/sigma
//   no_exit      same as current, but the cutoff is checked once after the loop over CVs
//   two_pass     dimension-major arrays; a branch-free pass computes the squared distance to every
//                kernel, a second pass evaluates exp and derivatives for kernels within the cutoff
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

struct Value {  // the relevant part of PLMD::Value
  enum {unset,periodic,notperiodic} periodicity=periodic;
  double max_minus_min=2*M_PI, inv_max_minus_min=1/(2*M_PI);
  static double pbc(double x) {  // Tools::pbc, round_toward_zero branch
    x+=100.0;
    return x>=0 ? x-int(x+0.5) : x-int(x-0.5);
  }
  double difference(double d1,double d2) const {
    if(periodicity==notperiodic) return d2-d1;
    double s=(d2-d1)*inv_max_minus_min;
    s=pbc(s);
    return s*max_minus_min;
  }
};

struct Kernel {
  double height;
  std::vector<double> center, sigma, inv_sigma;
};

struct Opes {
  unsigned ncv=0;
  double cutoff2=0, val_at_cutoff=0;
  std::vector<Value*> arguments;
  std::vector<Kernel> kernels;
  std::vector<double> period, inv_period;  // per CV, 0 if not periodic
  std::vector<double> dm_center, dm_inv_sigma, dm_height;  // dimension-major: dm_center[i*N+k]
  mutable std::vector<double> norm2_buf;
  mutable long n_within=0;

  double difference(unsigned i,double d1,double d2) const {return arguments[i]->difference(d1,d2);}

  double current(const std::vector<double>& x,std::vector<double>& der) const {
    std::vector<double> dist(ncv);
    double prob=0;
    for(const Kernel& G : kernels) {
      double norm2=0;
      bool within=true;
      for(unsigned i=0; i<ncv; i++) {
        dist[i]=difference(i,G.center[i],x[i])/G.sigma[i];
        norm2+=dist[i]*dist[i];
        if(norm2>=cutoff2) {within=false; break;}
      }
      if(!within) continue;
      const double der_val=G.height*std::exp(-0.5*norm2);
      prob+=der_val-G.height*val_at_cutoff;
      for(unsigned i=0; i<ncv; i++) der[i]-=dist[i]/G.sigma[i]*der_val;
    }
    return prob;
  }

  double inv_sigma(const std::vector<double>& x,std::vector<double>& der) const {
    std::vector<double> dist(ncv);
    double prob=0;
    for(const Kernel& G : kernels) {
      double norm2=0;
      bool within=true;
      for(unsigned i=0; i<ncv; i++) {
        dist[i]=difference(i,G.center[i],x[i])*G.inv_sigma[i];
        norm2+=dist[i]*dist[i];
        if(norm2>=cutoff2) {within=false; break;}
      }
      if(!within) continue;
      const double der_val=G.height*std::exp(-0.5*norm2);
      prob+=der_val-G.height*val_at_cutoff;
      for(unsigned i=0; i<ncv; i++) der[i]-=dist[i]*G.inv_sigma[i]*der_val;
    }
    return prob;
  }

  double no_exit(const std::vector<double>& x,std::vector<double>& der) const {
    std::vector<double> dist(ncv);
    double prob=0;
    for(const Kernel& G : kernels) {
      double norm2=0;
      for(unsigned i=0; i<ncv; i++) {
        dist[i]=difference(i,G.center[i],x[i])/G.sigma[i];
        norm2+=dist[i]*dist[i];
      }
      if(norm2>=cutoff2) continue;
      const double der_val=G.height*std::exp(-0.5*norm2);
      prob+=der_val-G.height*val_at_cutoff;
      for(unsigned i=0; i<ncv; i++) der[i]-=dist[i]/G.sigma[i]*der_val;
    }
    return prob;
  }

  double two_pass(const std::vector<double>& x,std::vector<double>& der) const {
    const unsigned N=dm_height.size();
    double* __restrict n2=norm2_buf.data();
    std::fill(n2,n2+N,0.);
    for(unsigned i=0; i<ncv; i++) {
      const double* __restrict c=&dm_center[i*N];
      const double* __restrict is=&dm_inv_sigma[i*N];
      const double xi=x[i], P=period[i], iP=inv_period[i];
      if(P>0) {
        for(unsigned k=0; k<N; k++) {
          double d=xi-c[k];
          d-=P*std::nearbyint(d*iP);
          const double t=d*is[k];
          n2[k]+=t*t;
        }
      } else {
        for(unsigned k=0; k<N; k++) {
          const double t=(xi-c[k])*is[k];
          n2[k]+=t*t;
        }
      }
    }
    double prob=0;
    for(unsigned k=0; k<N; k++) {
      if(n2[k]>=cutoff2) continue;
      n_within++;
      const double der_val=dm_height[k]*std::exp(-0.5*n2[k]);
      prob+=der_val-dm_height[k]*val_at_cutoff;
      for(unsigned i=0; i<ncv; i++) {
        double d=x[i]-dm_center[i*N+k];
        if(period[i]>0) d-=period[i]*std::nearbyint(d*inv_period[i]);
        const double is=dm_inv_sigma[i*N+k];
        der[i]-=d*is*is*der_val;
      }
    }
    return prob;
  }
};

static std::vector<std::vector<double>> read_rows(const char* fn) {
  std::ifstream f(fn);
  std::vector<std::vector<double>> rows;
  std::string line;
  while(std::getline(f,line)) {
    if(line.empty() || line[0]=='#') continue;
    std::istringstream s(line);
    std::vector<double> r;
    double v;
    while(s>>v) r.push_back(v);
    rows.push_back(r);
  }
  return rows;
}

int main(int argc,char** argv) {
  if(argc<3) {
    std::fprintf(stderr,"usage: %s STATE COLVAR [repeats]\n",argv[0]);
    return 1;
  }
  const auto state=read_rows(argv[1]);
  const auto colvar=read_rows(argv[2]);
  const int repeats=argc>3?std::atoi(argv[3]):3;
  Opes o;
  o.ncv=(state[0].size()-2)/2;  // time, centers, sigmas, height
  std::ifstream f(argv[1]);
  for(std::string line; std::getline(f,line);)
    if(line.find("SET kernel_cutoff")!=std::string::npos) {
      const double rc=std::stod(line.substr(line.rfind(' ')));
      o.cutoff2=rc*rc;
    }
  o.val_at_cutoff=std::exp(-0.5*o.cutoff2);
  static Value torsion;
  o.arguments.assign(o.ncv,&torsion);
  o.period.assign(o.ncv,2*M_PI);
  o.inv_period.assign(o.ncv,1/(2*M_PI));
  const unsigned N=state.size();
  o.dm_center.resize(N*o.ncv);
  o.dm_inv_sigma.resize(N*o.ncv);
  o.dm_height.resize(N);
  o.norm2_buf.resize(N);
  for(unsigned k=0; k<N; k++) {
    const auto& r=state[k];
    Kernel G;
    G.height=r[1+2*o.ncv];
    for(unsigned i=0; i<o.ncv; i++) {
      G.center.push_back(r[1+i]);
      G.sigma.push_back(r[1+o.ncv+i]);
      G.inv_sigma.push_back(1/r[1+o.ncv+i]);
      o.dm_center[i*N+k]=G.center[i];
      o.dm_inv_sigma[i*N+k]=G.inv_sigma[i];
    }
    o.dm_height[k]=G.height;
    o.kernels.push_back(G);
  }
  std::vector<std::vector<double>> points;
  for(const auto& r : colvar) points.emplace_back(r.begin()+1,r.begin()+1+o.ncv);
  std::printf("%u kernels, %u CVs, %zu points, best of %d\n",N,o.ncv,points.size(),repeats);

  using method=double (Opes::*)(const std::vector<double>&,std::vector<double>&) const;
  const std::vector<std::pair<const char*,method>> variants= {
    {"current",&Opes::current},{"inv_sigma",&Opes::inv_sigma},{"no_exit",&Opes::no_exit},{"two_pass",&Opes::two_pass}
  };
  double t_ref=0, p_ref=0, d_ref=0;
  for(const auto& v : variants) {
    double best=std::numeric_limits<double>::max(), sum_p=0, sum_d=0;
    o.n_within=0;
    for(int rep=0; rep<repeats; rep++) {
      sum_p=0;
      sum_d=0;
      const auto t0=std::chrono::steady_clock::now();
      for(const auto& x : points) {
        std::vector<double> der(o.ncv,0.);
        sum_p+=(o.*v.second)(x,der);
        for(double d : der) sum_d+=std::abs(d);
      }
      best=std::min(best,std::chrono::duration<double>(std::chrono::steady_clock::now()-t0).count());
    }
    if(t_ref==0) {
      t_ref=best;
      p_ref=sum_p;
      d_ref=sum_d;
    }
    const double rel=std::max(std::abs(sum_p/p_ref-1),std::abs(sum_d/d_ref-1));
    std::printf("%-10s %7.3f s  %6.2f ns/kernel  speedup %4.2f  relative change of the sums %.1e\n",
                v.first,best,best/(double(points.size())*N)*1e9,t_ref/best,rel);
    if(o.n_within>0) std::printf("           kernels within the cutoff: %.1f%%\n",100.*o.n_within/repeats/(double(points.size())*N));
  }
}
