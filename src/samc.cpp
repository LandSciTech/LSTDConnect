#include <Rcpp.h>
#include <samc.h>

#include <cstddef>
#include <vector>
#include <iostream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#define _P_FOCAL_OPENMP_ENADLED 1
#else
#define _P_FOCAL_OPENMP_ENADLED 0
#endif

namespace samc{

template<bool SYMMETRIC>
inline void construct_cache(
    cache& ca,
    const std::vector<kernel_point_t>& kernel,
    const Rcpp::NumericMatrix& resistance,
    const Rcpp::NumericMatrix& fidelity,
    const Rcpp::NumericMatrix& absorption){

  ca.kernel_size = kernel.size();
  ca.nrow = resistance.nrow();
  ca.ncol = resistance.ncol();

  ca.movement_rate.clear();
  ca.movement_rate.resize(ca.kernel_size*ca.nrow*ca.ncol, {0});

  ca.absorption.assign(absorption.begin(),absorption.end());

  std::ptrdiff_t max_offset = 0;
  std::ptrdiff_t min_offset = 0;
  for(const auto& k_point : kernel){
    const std::ptrdiff_t offset = k_point.y_off + k_point.x_off*ca.nrow;
    ca.kernel.push_back(-offset);
    max_offset = std::max(max_offset, offset);
    min_offset = std::min(min_offset, offset);
  }

  ca.left_extra_cols  = std::max((ca.nrow-1-min_offset)/ca.nrow, {0});
  ca.right_extra_cols = std::max((ca.nrow-1+max_offset)/ca.nrow, {0});
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(std::size_t x = 0; x<ca.ncol; x++){
    for(std::size_t y = 0; y<ca.nrow; y++){
      //Rcpp::Rcout << x << ", " << y << ", " << ca.death_rate[y+x*ca.nrow] << "\n";
      //bool printing = x==1 && y==1;

      double weighted_sum = 0;

      for(const auto& k_point : kernel){

        const std::size_t k_x = x+k_point.x_off;
        const std::size_t k_y = y+k_point.y_off;

        if((k_x >= 0 && k_x < ca.ncol) && (k_y >= 0 && k_y < ca.nrow)){
          weighted_sum += k_point.num/(resistance[k_y + k_x*ca.nrow]+(SYMMETRIC*resistance[y + x*ca.nrow]));
          //if(printing) Rcpp::Rcout << "a:" << x << ", " << y << ", " <<  k_x << ", " <<  k_y << ", " <<  resistance[k_y + k_x*ca.nrow] << ", " <<  k_point.num << "\n";
        }
      }
      double t_scalar = 0;
      double t_fidelity = fidelity[y+x*ca.nrow];
      const double t_absorption = absorption[y+x*ca.nrow];
      if(weighted_sum == 0){
        t_scalar = 0;
        t_fidelity = 1.0-t_absorption;
      }else{
        t_scalar = (1.0-(t_fidelity+t_absorption))/weighted_sum;
      }

      const double l_scalar     = t_scalar;
      const double l_fidelity   = t_fidelity;

      //if(printing) Rcpp::Rcout << "b:" << weighted_sum << ", " << l_scalar << ", " << l_fidelity << ", " << t_absorption << "\n";

      //Rcpp::Rcout << "(" << x << "," << y << "," << scalar << ")\n";
      //can be simd
      for(std::size_t k = 0; k < ca.kernel_size; k++){
        const auto& k_point = kernel[k];
        const std::size_t k_x = x+k_point.x_off;
        const std::size_t k_y = y+k_point.y_off;

        if((k_x >= 0 && k_x < ca.ncol) && (k_y >= 0 && k_y < ca.nrow)){
          ca.movement_rate[k+(k_y+k_x*ca.nrow)*ca.kernel_size] = (l_scalar * k_point.num)/(resistance[k_y + k_x*ca.nrow]+(SYMMETRIC*resistance[y + x*ca.nrow])) + (k_point.x_off == 0 && k_point.y_off == 0)*l_fidelity;

          //if(printing) Rcpp::Rcout << "c:" << k_x << ", " <<  k_y << ", " <<  ca.movement_rate[k+(k_y+k_x*ca.nrow)*ca.kernel_size] << "\n";
        }
        //if no, leave it 0, 0 is the default value anyway
      }

    }
  }
}

} /* namespace samc  */


// [[Rcpp::export(cache_samc_cpp)]]
Rcpp::XPtr<samc::cache> cache_samc(
    const Rcpp::NumericMatrix& kernel,
    const Rcpp::NumericMatrix& resistance,
    const Rcpp::NumericMatrix& fidelity,
    const Rcpp::NumericMatrix& absorption,
    const bool symmetric){

  std::vector<samc::kernel_point_t> kv{};

  const std::ptrdiff_t k_nrow = kernel.nrow();
  const std::ptrdiff_t k_ncol = kernel.ncol();

  for(std::ptrdiff_t y=0; y<k_nrow; y++){
    for(std::ptrdiff_t x=0; x<k_ncol; x++){
      if(kernel[y+x*k_nrow] || (y-k_nrow/2 == 0 && x-k_ncol/2 == 0)){
        kv.push_back({y-k_nrow/2, x-k_ncol/2, kernel[y+x*k_nrow]});
        //Rcpp::Rcout << kernel[y+x*k_nrow] << "\n";
      }
    }
  }

  if(kv.size() == 0){kv.push_back({0,0,0.0});}

  //for(size_t i = 0; i<kv.size(); i++){
  //  Rcpp::Rcout << i << ", " << kv[i].x_off << ", " << kv[i].y_off << ", " << kv[i].num << "\n";
  //}

  samc::cache* ca = new samc::cache;
  if(symmetric){
    samc::construct_cache<true>(*ca, kv, resistance, fidelity, absorption);
  }else{
    samc::construct_cache<false>(*ca, kv, resistance, fidelity, absorption);
  }
  Rcpp::XPtr<samc::cache> xp(ca, true);
  /*
  Rcpp::Rcout << "{"
  for(size_t x = 0; x<ca->ncol; x++){
    for(size_t y = 0; y<ca->nrow; y++){
      size_t i = y+x*ca->nrow;
      Rcpp::Rcout << "{\"x\": " << x << ", \"y\": " << y << ", \"absorption\": " << ca->absorption[i] << ", \"movement_rate\": [";
      bool first = true;
      for(size_t j = 0; j<ca->kernel_size; j++){
        if(first){first = false;}else{
          Rcpp::Rcout << ", ";
        }
        Rcpp::Rcout << ca->movement_rate[j+i*ca->kernel_size];
      }
      Rcpp::Rcout << "]},\n";
    }
  }*/
  //Rcpp::Rcout << *ca << '\n';
  return xp;
}

inline void samc_one_step(
    const samc::cache& ca,
    const double* const pop_in,
    const double* const dead_in,
    double* const pop_out,
    double* const dead_out){

  //const auto indexes     = std::ranges::iota_view<size_t, size_t>(0, ca.death_rate.size());
  //const auto connections = std::ranges::iota_view<size_t, size_t>(0, ca.kernel_size);

  //const std::size_t offset = ca.nrow * ca.left_extra_cols;

  const double* const p_in = pop_in;// + (ca.nrow * ca.left_extra_cols);
  const double* const d_in = dead_in;

  double* const p_out = pop_out;// + (ca.nrow * ca.left_extra_cols);
  double* const d_out = dead_out;

  //Rcpp::Rcout << ca.death_rate.size() << " " << ca.kernel_size << " " << offset <<"\n";
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for(std::size_t i = 0; i<ca.absorption.size(); i++){
    d_out[i] = d_in[i]+ca.absorption[i]*p_in[i];
    double acc = 0;
    for(std::size_t con = 0; con < ca.kernel_size; con++){
      //Rcpp::Rcout << i << ", " << i%ca.nrow << ", " << i/ca.nrow << ", " << con << ", " << i*ca.kernel_size+con << ", " << i+ca.kernel[con] << ", " << ca.movement_rate[i*ca.kernel_size+con] << '\n';
      acc += ca.movement_rate[i*ca.kernel_size+con]*p_in[i+ca.kernel[con]];
    }
    p_out[i] = acc;
  }
}
/*
// [[Rcpp::export(samc_print_cache_as_matrix)]]
void samc_print_cache_as_matrix(const Rcpp::XPtr<samc::cache>& ca){
  const samc::cache& c = *ca;
  Rcpp::Rcout << "[";
  for(size_t i = 0; i<c.nrow; i++){
    Rcpp::Rcout << "[";
    for(size_t j = 0; j<c.ncol; j++){
        //sorry, I ran out of time, good luck
    }
  }
}*/

// [[Rcpp::export(samc_cache_sizes_cpp)]]
std::vector<size_t> samc_cache_sizes(const Rcpp::XPtr<samc::cache>& ca){
  return {ca->nrow, ca->ncol, ca->left_extra_cols, ca->right_extra_cols};
}

// [[Rcpp::export(samc_step_cpp)]]
Rcpp::List samc_step(
    std::vector<long> steps,
    const Rcpp::XPtr<samc::cache>& ca,
    const Rcpp::NumericMatrix& pop_in,
    Rcpp::NumericMatrix dead_in){

  if(steps.size() <= 0){
    Rcpp::Rcerr << "We need at least one step number\n";
    return {};
  }

  std::sort(std::begin(steps), std::end(steps));
  steps.erase(std::unique(std::begin(steps), std::end(steps)), std::end(steps));
  //steps is now a serted list of step numbers with no duplicates

  if(steps[0] < 0){
    Rcpp::Rcerr << "We cannot step backwards, all step numbers must be at least 0\n";
    return {};
  }

  std::vector<double> pop_a(ca->nrow*(ca->ncol+ca->left_extra_cols+ca->right_extra_cols), 0.0);
  std::vector<double> pop_b(ca->nrow*(ca->ncol+ca->left_extra_cols+ca->right_extra_cols), 0.0);

  std::memcpy(&pop_a[ca->nrow*ca->left_extra_cols], &pop_in[0], ca->nrow*ca->ncol*sizeof(double));

  std::vector<double> dead_b(ca->nrow*ca->ncol, 0.0);

  std::vector<Rcpp::NumericMatrix> pops{};
  std::vector<Rcpp::NumericMatrix> deads{};

  double* const p_a = &pop_a[ca->nrow*ca->left_extra_cols];
  double* const p_b = &pop_b[ca->nrow*ca->left_extra_cols];
  double* const d_a = &dead_in[0];
  double* const d_b = &dead_b[0];

  const double* p_in = p_a;
  const double* d_in = d_a;
  double* p_out = p_b;
  double* d_out = d_b;

  long last_i = 0;
  for(long i : steps){
    for(long j=0; j<i-last_i; j++){
      //step once
      samc_one_step(*ca, p_in, d_in, p_out, d_out);

      p_in = p_out;
      d_in = d_out;
      if(p_out == p_a){
        p_out = p_b;
        d_out = d_b;
      }else{
        p_out = p_a;
        d_out = d_a;
      }
    }
    pops.emplace_back(int(ca->nrow), int(ca->ncol));
    deads.emplace_back(int(ca->nrow), int(ca->ncol));

    std::memcpy(&pops.back()[0], p_in, ca->nrow*ca->ncol*sizeof(double));
    std::memcpy(&deads.back()[0], d_in, ca->nrow*ca->ncol*sizeof(double));

    //add to output lists
  }

  return Rcpp::List::create(Rcpp::Named("time") = steps, Rcpp::Named("occ") = pops, Rcpp::Named("deaths") = deads);
}







