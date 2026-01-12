/* Benchmark of boost::hub against plf::hive.
 * 
 * Copyright 2026 Joaquin M Lopez Munoz.
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 */

#include <algorithm>
#include <array>
#include <chrono>
#include <numeric>

std::chrono::high_resolution_clock::time_point measure_start, measure_pause;

template<typename F>
double measure(F f)
{
  using namespace std::chrono;

  static const int              num_trials = 10;
  static const milliseconds     min_time_per_trial(200);
  std::array<double,num_trials> trials;

  for(int i = 0; i < num_trials; ++i) {
    int                               runs = 0;
    high_resolution_clock::time_point t2;
    volatile decltype(f())            res; /* to avoid optimizing f() away */

    measure_start = high_resolution_clock::now();
    do{
      res = f();
      ++runs;
      t2 = high_resolution_clock::now();
    }while(t2 - measure_start<min_time_per_trial);
    trials[i] =
      duration_cast<duration<double>>(t2 - measure_start).count() / runs;
  }

  std::sort(trials.begin(), trials.end());
  return std::accumulate(
    trials.begin() + 2, trials.end() - 2, 0.0)/(trials.size() - 4);
}

void pause_timing()
{
  measure_pause = std::chrono::high_resolution_clock::now();
}

void resume_timing()
{
  measure_start += std::chrono::high_resolution_clock::now() - measure_pause;
}

#include <boost/core/detail/splitmix64.hpp>
#include <boost/hub.hpp>
#include <cmath>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <plf_hive.h>
#include <stdexcept>
#include <string>
#include <vector>

struct element
{
  element(int n_): n{n_} {}

#if defined(NONTRIVIAL_ELEMENT)
  element(element&& x): n{x.n}
  {
    std::memcpy(payload, x.payload, sizeof(payload));
    std::memset(x.payload, 0, sizeof(payload));
  }

  element& operator=(element&& x)
  {
    n = x.n;
    std::memcpy(payload, x.payload, sizeof(payload));
    std::memset(x.payload, 0, sizeof(payload));
    return *this;
  }
#endif

  operator int() const { return n; }

  int n;
  std::string str = std::to_string(n);
  char payload[ELEMENT_SIZE - sizeof(int)];
};

struct urbg
{
  using result_type = boost::uint64_t;

  static constexpr result_type min() { return 0; }
  static constexpr result_type max() { return (result_type)(-1); }

  urbg() = default;
  explicit urbg(result_type seed): rng{seed} {}

  result_type operator()() { return rng(); }

  boost::detail::splitmix64 rng;
};

template<typename Container, typename Iterator>
void erase_void(Container& x, Iterator it)
{
  x.erase(it);
}

template<typename... Args, typename Iterator>
void erase_void(boost::hub<Args...>& x, Iterator it)
{
  x.erase_void(it);
}

template<typename Container>
Container make(std::size_t n, double erasure_rate)
{
  std::size_t   m = (std::size_t)((double)n / (1.0 - erasure_rate));
  std::uint64_t erasure_cut = 
    (std::uint64_t)(erasure_rate * (double)(std::uint64_t)(-1));
  std::size_t   reinsertion_count =
    (std::size_t)(erasure_rate * 0.5 * (double)n);

  Container                                 c;
  urbg                                      rng;
  std::vector<typename Container::iterator> iterators;

  iterators.reserve(m);
  for(std::size_t i = 0; i < m; ++i) iterators.push_back(c.insert((int)rng()));
  std::shuffle(iterators.begin(), iterators.end(), rng);
  for(auto it: iterators) {
    if(rng() < erasure_cut) erase_void(c, it);
  }
  for(std::size_t i = 0; i < reinsertion_count; ++i) c.insert((int)rng());
  return c;
}

template<typename FHive, typename FHub>
void benchmark(const char* title, FHive fhive, FHub fhub)
{
  static constexpr std::size_t size_limit =
    sizeof(std::size_t) == 4?  800ull * 1024ull * 1024ull:
                              2048ull * 1024ull * 1024ull;

  std::cout << std::string(41, '-') << "\n"
            << title << "\n"
            << "sizeof(element): " << sizeof(element) << "\n";
  std::cout << std::left << std::setw(11) << "" << "container size\n" << std::right
            << std::left << std::setw(11) << "erase rate" << std::right;
  for(std::size_t i = 3; i <= 7; ++i)
  {
    std::cout << "1.E" << i << " ";
  }
  std::cout << std::endl;
  for(double erasure_rate = 0.0; erasure_rate <= 0.9; erasure_rate += 0.1) {
    std::cout << std::left << std::setw(11) << erasure_rate << std::right << std::flush;
    for(std::size_t i = 3; i <= 7; ++i) {
      std::size_t n = (std::size_t)std::pow(10.0, (double)i);
      if((double)n * (double)sizeof(element) / (1.0 - erasure_rate) > (double)size_limit) {
        std::cout << "---- " << std::flush;
        continue;
      }

      auto thive = measure([&] { return fhive(n, erasure_rate); });
      auto thub = measure([&] { return fhub(n, erasure_rate); });
      std::cout << std::fixed << std::setprecision(2)
                << thive / thub << " "
                << std::defaultfloat << std::setprecision(6) << std::flush;
    }
    std::cout << std::endl;
  }
}

template<typename Container>
struct create
{
  unsigned int operator()(std::size_t n, double erasure_rate) const
  {
    unsigned int res = 0;
    {
      auto c = make<Container>(n, erasure_rate);
      res = c.size();
      pause_timing();
    }
    resume_timing();
    return res;
  }
};

template<typename Container>
struct create_and_destroy
{
  unsigned int operator()(std::size_t n, double erasure_rate) const
  {
    auto c = make<Container>(n, erasure_rate);
    return c.size();
  }
};

template<typename Container>
struct prepare
{
  const Container& get_container(std::size_t n_, double erasure_rate_)
  {
    if(n_ != n || erasure_rate_ != erasure_rate) {
      pause_timing();
      n = n_;
      erasure_rate = erasure_rate_;
      c.clear();
      //c.shrink_to_fit();
      c = make<Container>(n, erasure_rate);
      resume_timing();
    }
    return c;
  }

  std::size_t n = 0;
  double      erasure_rate = 0.0;
  Container   c;
};

template<typename Container>
struct for_each: prepare<Container>
{
  unsigned int operator()(std::size_t n, double erasure_rate)
  {
    unsigned int res = 0;
    auto& c = this->get_container(n, erasure_rate);
    for(const auto& x: c) res += (unsigned int)x; 
    return res;
  }
};

template<typename Container>
struct visit_all: prepare<Container>
{
  unsigned int operator()(std::size_t n, double erasure_rate)
  {
    unsigned int res = 0;
    auto& c = this->get_container(n, erasure_rate);
    c.visit_all([&] (const auto& x) { res += (unsigned int)x; });
    return res;
  }
};

int main()
{
  try{
    using hive = plf::hive<element>;
    using hub = boost::hub<element>;

    benchmark(
      "insert, erase, insert", 
      create<hive>{}, create<hub>{});
    benchmark(
      "insert, erase, insert, destroy", 
      create_and_destroy<hive>{}, create_and_destroy<hub>{});
    benchmark(
      "for_each", 
      for_each<hive>{}, for_each<hub>{});
    benchmark(
      "visit_all", 
      for_each<hive>{}, visit_all<hub>{});
  }
  catch(const std::exception& e) {
    std::cerr << e.what() << std::endl;
  }
}