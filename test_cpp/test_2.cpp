#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>

double EPS = 0.01;
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





	double golden_ratio(std::function<double(double)> f, const double& a, const double& b, double eps) {
		double a_k = a;
		double b_k = b;
		size_t k = 1;
		double x_k = a_k + (1 - TAU) * (b_k - a_k);
		double y_k = a_k + (b_k - a_k) * TAU;
		double f_x = f(x_k);
		double f_y = f(y_k);
		while (b_k - a_k > eps) {
			if (f_x > f_y) {
				//step 2
				a_k = x_k; //k++ = k
				x_k = y_k;
				f_x = f_y;
				y_k = a_k + TAU * (b_k - a_k);
				f_y = f(y_k);
			}
			else {
				//step 3
				a_k = a_k; // k++ = k
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

		size_t k = 0;
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
		std::cout << golden_ratio(f,-1, 1, 0.01) << std::endl;
	}


	void test_fib() {
		std::cout << method_fib(-1, 1, 0.01) << std::endl;
	}
}





namespace lab_3 {
	struct Point {
		double x1;
		double x2;

		Point(double a, double b) : x1(a), x2(b) {};

		Point operator+(const Point& other) const {
			return { x1 + other.x1, x2 + other.x2 };
		}

		Point operator-(const Point& other) const {
			return { x1 - other.x1, x2 - other.x2 };
		}

		Point operator-() const {
			return { -x1, -x2 };
		}

		Point operator*(double num) const {
			return { num * x1, num * x2 };
		}

		friend Point operator*(double scalar, const Point& p) {
			return { scalar * p.x1, scalar * p.x2 };
		}

		double dot(const Point& other) const {
			return x1 * other.x1 + x2 * other.x2;
		}

		double product() const {
			return x1 * x1 + x2 * x2;
		}

		double norm() const {
			return std::sqrt(product());
		}
	};

	double f(const Point& p) {
		return std::exp(3 * p.x1) + std::pow(p.x1 + p.x2, 2) + std::exp(2* p.x2);
	}

	Point grad_f(const Point& p) {
		return { 3 * std::exp(3 * p.x1) + 2 * (p.x1 + p.x2), 2 * (p.x1 + p.x2) + 2 * std::exp(2 * p.x2) };
	}

	Point grad_armijo(Point& p, double (*f)(const Point&), Point(*grad_f)(const Point&),
		double alpha, double gamma, double theta, double eps)
	{
		Point x_k = p;
		long long k = 0;

		while (grad_f(x_k).norm() > eps) {
			k += 1;
			double alpha_k = alpha;

			Point grad = grad_f(x_k);
			while (f(x_k - alpha_k * grad) > f(x_k) - gamma * alpha_k * grad.dot(grad)) {
				alpha_k *= theta;
				
			}
			
			x_k = x_k - alpha_k * grad;
		}
		std::cout << 'K = ' << k << ' ';
		return x_k;
	}

	void test_dif_params() {
		Point x0(1.0, 1.0);

		std::cout << "Standart TEST" << std::endl;
		auto result = grad_armijo(x0, f, grad_f, 1, 0.5, 0.5, 0.01);
		std::cout << result.x1 << " " << result.x2 << std::endl;
		
		
		std::cout << "changing alpha" << std::endl;
		std::vector<double> alphas = { 1.0, 0.5, 0.1, 0.01 };
		for (auto kent : alphas) {
			auto result = grad_armijo(x0, f, grad_f, kent, 0.5, 0.5, 0.01);
			std::cout << "#" << kent << "# " << result.x1 << " " << result.x2 << std::endl;
		}

		std::cout << "changing gamma" << std::endl;
		std::vector<double> gammas = { 0.1, 0.3, 0.5, 0.7 };
		for (auto kent : gammas) {
			auto result = grad_armijo(x0, f, grad_f, 1, kent, 0.5, 0.01);
			std::cout << "#" << kent << "# " << result.x1 << " " << result.x2 << std::endl;
		}


		std::cout << "changing thrta" << std::endl;
		std::vector<double> thetas = { 0.1, 0.3, 0.4, 0.5, 0.8 };
		for (auto kent : thetas) {
			auto result = grad_armijo(x0, f, grad_f, 1, 0.5, kent, 0.01);
			std::cout << "#" << kent << "# " << result.x1 << " " << result.x2 << std::endl;
		}

		
		


	}
}


namespace lab_4 {
	using namespace lab_3;
	using namespace lab_2;

	Point grad_faster(Point& x0, double (*f)(const Point&), Point(*grad_f)(const Point&),
		double eps, double alpha_max = 10.0)
	{

		Point x_k = x0;
		int k = 0;

		while (grad_f(x_k).norm() > eps) {
			k++;

			Point grad = grad_f(x_k);

			auto F_alpha = [&grad, &x_k, &f](double alpha) -> double {
				Point z = x_k - alpha * grad;
				return f(z);
				};

			double alpha_k = golden_ratio(F_alpha, 0, alpha_max, eps);

			x_k = x_k - alpha_k * grad;
		}
		std::cout << 'K = ' << k << ' ';
		return x_k;
	}


	void test() {
		lab_3::Point x0(1.0, 1.0);
		std::cout << "GRAD_FAStER" << std::endl;
		auto result = lab_4::grad_faster(x0, lab_3::f, lab_3::grad_f, 0.01);
		std::cout << result.x1 << " " << result.x2  << std::endl;
	}
}



int main() {
	//std::cout << TAU;
	//lab_2::test_golden();
	//lab_2::test_fib();
	lab_3::test_dif_params();
	std::cout << std::endl;
	lab_4::test();

	return 0;
}
