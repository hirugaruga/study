#include <iostream>
#include <vector>
#include <algorithm>

std::vector<int> sequentialDigits(int low, int high) {
  std::string c = "123456789";
  std::vector<int> a;

  for (int i = 0; i < c.size(); i++) {
	for (int j = i + 1; j <= c.size(); j++) {
	  int curr = std::stoi(c.substr(i, j - i));
	  if (curr <= high && curr >= low) {
		a.push_back(curr);
	  }
	}
  }
  std::sort(a.begin(), a.end());
  return a;
}

  int main() {
	int low;
	int high;
	std::cin >> low >> high;
	std::vector<int> ans = sequentialDigits(low, high);
	for (auto x : ans) {
	  std::cout << x << ' ';
	}
  }

#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>


double EPS = 0.001;
double TAU = (std::sqrt(5) - 1) / 2;

namespace lab_2 {
  double f(double x) {
    return std::exp(x * x + x) + std::exp(x * x + 1);
  }

  std::vector<double> fibon(int n) {
    std::vector<double> fib;
    if (n <= 0) return fib;
    fib.reserve(n);
    if (n >= 1) fib.push_back(1);
    if (n >= 2) fib.push_back(1);

    for (size_t i = 2; i < n; i++) {
      fib.push_back(fib[i - 1] + fib[i - 2]);
    }

    return fib;
  }





  double golden_ratio(const double& a, const double& b, double eps) {
    double a_k = a;
    double b_k = b;
    size_t k = 1;
    double x_k = a_k + (1 - TAU) * (b_k - a_k);
    double y_k = a_k + TAU * (b_k - a_k);
    double f_x = f(x_k);
    double f_y = f(y_k);
    while (b_k - a_k > eps) {
      if (f_x > f_y) {
        //step 2
        a_k = x_k; //k++ = k
        b_k = b_k;
        x_k = y_k;
        f_x = f_y;
        double y_k = a_k + TAU * (b_k - a_k);
        f_y = f(y_k);
      }
      else {
        //step 3
        a_k = x_k; // k++ = k
        b_k = y_k;
        y_k = x_k;
        f_y = f_x;
        x_k =  a_k + (1 - TAU) * (b_k - a_k);
        f_x = f(x_k);
       }
      k++;
    }
    return (a_k + b_k) / 2;

  }

  long long n_founder(const double& a, const double& b, std::vector<double>& fibs, const double& eps) {
    long long k = 1;
    while ((b - a) / fibs[k] > eps) k++;
    return k;
  }


  double method_fib(const double& a, const double& b, double eps) {
    std::vector<double> F = fibon(1000);
    double alpha = eps * 0.1;
    double a_k = a;
    double b_k = b;
    size_t k = 1;
    size_t n = n_founder(a, b, F, eps);
    double x_k = a_k + (F[n - k - 1] / F[n - k + 1]) * (b_k - a_k);
    double y_k = a_k + (F[n - k] / F[n - k + 1]) * (b_k - a_k);
    double f_x = f(x_k);
    double f_y = f(y_k);
    while (k < n - 1){
      if (f_x > f_y) {
        a_k = x_k;
        //b_k = b_k;
        x_k = y_k;
        f_x = f_y;
        y_k = a_k + (F[n - k] / F[n - k + 1]) * (b_k - a_k);
        f_y = f(y_k);
      }
      else {
        //a_k = a_k;
        b_k = y_k;
        y_k = x_k;
        f_y = f_x;
        x_k = a_k + (F[n - k - 1] / F[n - k + 1]) * (b_k - a_k);
        f_x = f(x_k);
      }
      k++;

      double x_n = x_k;
      double y_n = y_k + alpha;

      double f_x_n = f(x_n);
      double f_y_n = f(y_n);

      if (f_x_n > f_y_n) {
        a_k = x_n;
      }
      else {
        b_k = y_n;
      }
      return (a_k + b_k) / 2;
    } 

  }




  void test_golden() {
    std::cout << golden_ratio(-1, 1, 0.01) << std::endl;
  }


  void test_fib() {
    std::cout << method_fib(-1, 1, 0.01) << std::endl;
  }
}


int main() {
  lab_2::test_golden();
  lab_2::test_fib();

  return 0;
}