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

inline void construct_cache(cache& ca, const std::vector<kernel_point_t>& kernel, const Rcpp::NumericMatrix& permiability, const Rcpp::NumericMatrix& death_rate){
  
  ca.kernel_size = kernel.size();
  ca.nrow = permiability.nrow();
  ca.ncol = permiability.ncol();
  
  if((ca.ncol * ca.nrow)==0){
    throw "There are no cells, nrow*ncol == 0";
  }
  
  if((std::size_t)(death_rate.ncol()) != ca.ncol || ((std::size_t)death_rate.nrow()) != ca.nrow){
    throw "death_rate's size does not match permiability's size";
  }
  
  ca.movement_rate.clear();
  ca.movement_rate.resize(ca.kernel_size*ca.nrow*ca.ncol, {0});
  ca.death_rate.assign(death_rate.begin(),death_rate.end());
  
  
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
  
  #pragma omp parallel for 
  for(std::size_t x = 0; x<ca.ncol; x++){
    for(std::size_t y = 0; y<ca.nrow; y++){
      //std::cout << x << ", " << y << ", " << ca.death_rate[y+x*ca.nrow] << "\n";
      double weighted_sum = 0;
      
      for(const auto& k_point : kernel){
        //const auto& k_point = kernel[k];
        const std::size_t k_x = x+k_point.x_off;
        const std::size_t k_y = y+k_point.y_off;
        //std::cout << "(" << x << ", " << y << ", " << k << ", ";
        if((k_x >= 0 && k_x < ca.ncol) && (k_y >= 0 && k_y < ca.nrow)){
          //std::cout << permiability[k_y + k_x*ca.nrow] * k_point.num <<")\n";
          weighted_sum += permiability[k_y + k_x*ca.nrow] * k_point.num;
        }//else{
        //std::cout << "0)\n";
        // weighted_sum += 0.0;
        //}
      }
      
      const double scalar = (weighted_sum)?((1.0-death_rate[y+x*ca.nrow])/weighted_sum):(0);
      //std::cout << "(" << x << "," << y << "," << scalar << ")\n";
      //can be simd
      for(std::size_t k = 0; k < ca.kernel_size; k++){
        const auto& k_point = kernel[k];
        const std::size_t k_x = x+k_point.x_off;
        const std::size_t k_y = y+k_point.y_off;
        
        if((k_x >= 0 && k_x < ca.ncol) && (k_y >= 0 && k_y < ca.nrow)){
          ca.movement_rate[k+(k_y+k_x*ca.nrow)*ca.kernel_size] = scalar * permiability[k_y + k_x*ca.nrow];
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
    const Rcpp::NumericMatrix& permiability, 
    const Rcpp::NumericMatrix& death_rate) {
  
  std::vector<samc::kernel_point_t> kv{};
  
  const std::ptrdiff_t k_nrow = kernel.nrow();
  const std::ptrdiff_t k_ncol = kernel.ncol();
  for(std::ptrdiff_t y=0; y<k_nrow; y++){
    for(std::ptrdiff_t x=0; x<k_ncol; x++){
      if(kernel[y+x*k_nrow]){
        kv.push_back({y-k_nrow/2, x-k_ncol/2, kernel[y+x*k_nrow]});
      }
    }
  }
  
  samc::cache* ca = new samc::cache;
  samc::construct_cache(*ca, kv, permiability, death_rate);
  Rcpp::XPtr<samc::cache> xp(ca, true);
  return xp;
}

inline void samc_step(
    const samc::cache& ca,
    const Rcpp::NumericMatrix& pop_in,
    const Rcpp::NumericMatrix& dead_in,
    Rcpp::NumericMatrix& pop_out,
    Rcpp::NumericMatrix& dead_out){
  
  //const auto indexes     = std::ranges::iota_view<size_t, size_t>(0, ca.death_rate.size());
  //const auto connections = std::ranges::iota_view<size_t, size_t>(0, ca.kernel_size);
  
  const std::size_t offset = ca.nrow * ca.left_extra_cols;
  
  //std::cout << ca.death_rate.size() << " " << ca.kernel_size << " " << offset <<"\n";
  #pragma omp parallel for 
  for(std::size_t i = 0; i<ca.death_rate.size(); i++){
    dead_out[i+offset] = dead_in[i+offset]+ca.death_rate[i]*pop_in[i+offset];
    double acc = 0;
    for(std::size_t con = 0; con < ca.kernel_size; con++){
      //std::cout << i << "\t" << con << "\t" << i*ca.kernel_size+con << "\t" << i+offset+ca.kernel[con]  <<"\n";
      //return ca.movement_rate[i][con]*pop_in[i+offset+ca.kernel[con]];
      acc += ca.movement_rate[i*ca.kernel_size+con]*pop_in[i+offset+ca.kernel[con]];
      //return 1.0;
    }
  pop_out[i+offset] = acc;
  }
}

// [[Rcpp::export(samc_cache_sizes_cpp)]]
std::vector<size_t> samc_cache_sizes(const Rcpp::XPtr<samc::cache>& ca){
  return {ca->nrow, ca->ncol, ca->left_extra_cols, ca->right_extra_cols};
}

// [[Rcpp::export(samc_print_cache_cpp)]]
void samc_print_cache(const Rcpp::XPtr<samc::cache>& ca){
  std::cout << *ca;
}


//Rcpp::SubMatrix<REALSXP> samc_one_step(
// [[Rcpp::export(samc_step_cpp)]]
Rcpp::List samc_step(
    std::vector<long> steps,
    const Rcpp::XPtr<samc::cache>& ca,
    Rcpp::NumericMatrix& pop_in,
    Rcpp::NumericMatrix& dead_in){
  
  if(steps.size() <= 1){
    std::cerr << "We need at least one step number\n";
    return {};
  }
  
  std::sort(std::begin(steps), std::end(steps));
  steps.erase(std::unique(std::begin(steps), std::end(steps)), std::end(steps));
  //steps is now a serted list of step numbers with no duplicates
  
  if(steps[0] < 0){
    std::cerr << "We cannot step backwards, all step numbers must be at least 0\n";
    return {};
  }
  
  Rcpp::NumericMatrix pop_in_b( ca->nrow, (ca->ncol + ca->left_extra_cols + ca->right_extra_cols));
  Rcpp::NumericMatrix dead_in_b(ca->nrow, (ca->ncol + ca->left_extra_cols + ca->right_extra_cols));
  
  const Rcpp::NumericMatrix* last_pop_p = &pop_in;
  const Rcpp::NumericMatrix* last_dead_p = &dead_in;
  
  Rcpp::NumericMatrix* next_pop_p = &pop_in_b;
  Rcpp::NumericMatrix* next_dead_p = &dead_in_b;
  
  
  std::vector<Rcpp::NumericMatrix> pops{};
  std::vector<Rcpp::NumericMatrix> deads{};
  
  long last_i = 0;
  
  for(size_t index=0; index<steps.size(); index++){
    long i = steps[index];
    for(long j = 0; j<i-last_i; j++){
      samc_step(*ca, *last_pop_p, *last_dead_p, *next_pop_p, *next_dead_p);
      last_pop_p = next_pop_p;
      last_dead_p = next_dead_p;
      if(next_pop_p == &pop_in_b){
        next_pop_p = &pop_in;
        next_dead_p = &dead_in;
      }else{
        next_pop_p = &pop_in_b;
        next_dead_p = &pop_in_b;
      }
    }
    last_i = i;
    if(index+1 < steps.size()){
      pops.emplace_back(Rcpp::clone(*last_pop_p));
      deads.emplace_back(Rcpp::clone(*last_dead_p));
    }else{
      pops.emplace_back(*last_pop_p);
      deads.emplace_back(*last_dead_p);
    }
  }
  
  return Rcpp::List::create(Rcpp::Named("steps") = steps, Rcpp::Named("population") = pops, Rcpp::Named("deaths") = deads);
}







