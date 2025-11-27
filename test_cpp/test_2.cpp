#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include <functional>

double EPS = 0.01;
double TAU = (std::sqrt(5) - 1) / 2;

using namespace std;


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
				x_k = a_k + (1 - TAU) * (b_k - a_k);
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
		while (k < n - 1) {
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
		std::cout << golden_ratio(f, -1, 1, 0.01) << std::endl;
	}


	void test_fib() {
		std::cout << method_fib(-1, 1, 0.01) << std::endl;
	}
}

namespace lab_3 {
	struct Point {
		std::vector<double> coor;

		Point(const std::vector<double>& temp) : coor(temp) {};
		Point(int n, const std::vector<double>& temp) : coor(temp) {};
		Point(double a, double b) : coor({ a,b }) {};

		Point() = default;


		//Point operator-(const Point& other) const {
		//	return { x1 - other.x1, x2 - other.x2 };
		//}

		Point operator-(const Point& other) const {
			auto temp = coor;
			for (int i = 0; i < coor.size(); ++i) {
				temp[i] = coor[i] - other.coor[i];
			}

			return temp;
		}

		Point operator*(double num) const {
			auto temp = coor;
			for (size_t i = 0; i < coor.size(); i++) {
				temp[i] *= num;
			}
			return temp;



		}

		double operator[](size_t index) {
			return coor[index];
		}

		Point operator+(const Point& other) const {
			auto temp = coor;
			for (int i = 0; i < coor.size(); ++i) {
				temp[i] = coor[i] + other.coor[i];
			}

			return temp;
		}

		friend Point operator*(double scalar, const Point& p) {
			auto val = p.coor;
			for (size_t i = 0; i < p.coor.size(); ++i) {
				val[i] *= scalar;
			}
			return val;
		}

		double dot(const Point& other) const {
			//return x1 * other.x1 + x2 * other.x2;
			double val = 0;
			for (size_t i = 0; i < coor.size(); ++i) {
				val += coor[i] * other.coor[i];
			}
			return val;
		}

		double product() const {
			double val = 0;
			for (size_t i = 0; i < coor.size(); ++i) {
				val += coor[i] * coor[i];
			}
			return val;
		}

		double norm() const {
			return std::sqrt(product());
		}
	};

	double f(const Point& p) {
		return std::exp(3 * p.coor[0]) + std::pow(p.coor[0] + p.coor[1], 2) + std::exp(2 * p.coor[1]);
	}

	Point grad_f(const Point& p) {
		return { 3 * std::exp(3 * p.coor[0]) + 2 * (p.coor[0] + p.coor[1]) , 2 * (p.coor[0] + p.coor[1]) + 2 * std::exp(2 * p.coor[1]) };
	}

	Point grad_armijo(Point& p, double (*f)(const Point&), Point(*grad_f)(const Point&),
		double alpha, double gamma, double theta, double eps)
	{
		Point x_k = p;
		long long k = 0;
		Point z = grad_f(x_k);
		double zz = z.norm();
		while (grad_f(x_k).norm() > eps) {
			k += 1;
			double alpha_k = alpha;

			Point grad = grad_f(x_k);
			while (f(x_k - alpha_k * grad) > f(x_k) - gamma * alpha_k * grad.dot(grad)) {
				alpha_k *= theta;

			}

			x_k = x_k - alpha_k * grad;
		}
		std::cout << "K = " << k << ' ';
		return x_k;
	}

	void test_dif_params() {
		Point x0(1, 1);

		Point z = grad_f(x0);
		double zz = z.norm();
		std::cout << "Standart TEST" << std::endl;
		auto result = grad_armijo(x0, f, grad_f, 1, 0.5, 0.5, 0.01);
		std::cout << result.coor[0] << " " << result.coor[1] << std::endl;


		std::cout << "changing alpha" << std::endl;
		std::vector<double> alphas = { 1.0, 0.5, 0.1, 0.01 };
		for (auto kent : alphas) {
			auto result = grad_armijo(x0, f, grad_f, kent, 0.5, 0.5, 0.01);
			std::cout << "#" << kent << "# " << result.coor[0] << " " << result.coor[1] << std::endl;
		}

		std::cout << "changing gamma" << std::endl;
		std::vector<double> gammas = { 0.1, 0.3, 0.5, 0.7 };
		for (auto kent : gammas) {
			auto result = grad_armijo(x0, f, grad_f, 1, kent, 0.5, 0.01);
			std::cout << "#" << kent << "# " << result.coor[0] << " " << result.coor[1] << std::endl;
		}


		std::cout << "changing thrta" << std::endl;
		std::vector<double> thetas = { 0.1, 0.3, 0.4, 0.5, 0.8 };
		for (auto kent : thetas) {
			auto result = grad_armijo(x0, f, grad_f, 1, 0.5, kent, 0.01);
			std::cout << "#" << kent << "# " << result.coor[0] << " " << result.coor[1] << std::endl;
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
			/*
			auto F_alpha = [&](double alpha) -> double {
				Point z = x_k - alpha * grad;
				return f(z);
			};*/

			double alpha_k = golden_ratio([&](double alpha) -> double { Point z = x_k - alpha * grad; return f(z);}, 0, alpha_max, eps);

			x_k = x_k - alpha_k * grad;
		}
		std::cout << "K = " << k << ' ';
		return x_k;
	}


	void test() {
		lab_3::Point x0(1, 1);
		std::cout << "GRAD_FAStER" << std::endl;
		auto result = lab_4::grad_faster(x0, lab_3::f, lab_3::grad_f, 0.01);
		std::cout << result.coor[0] << " " << result.coor[1] << std::endl;
	}
}

namespace lab_5{

	using namespace lab_2;
    class Point {
    public:
        vector<double> coor;

        Point(const vector<double>& temp = {}) : coor(temp) {};
        Point(double a, double b, double c) : coor({ a, b, c }) {};

        Point operator-(const Point& other) const {
            vector<double> temp;
            for (size_t i = 0; i < coor.size(); ++i) {
                temp.push_back(coor[i] - other.coor[i]);
            }
            return Point(temp);
        }

        Point operator+(const Point& other) const {
            vector<double> temp;
            for (size_t i = 0; i < coor.size(); ++i) {
                temp.push_back(coor[i] + other.coor[i]);
            }
            return Point(temp);
        }

        Point operator*(double num) const {
            vector<double> temp;
            for (size_t i = 0; i < coor.size(); ++i) {
                temp.push_back(coor[i] * num);
            }
            return Point(temp);
        }

        friend Point operator*(double scalar, const Point& p) {
            return p * scalar;
        }

        double operator[](size_t index) const {
            return coor[index];
        }

        double dot(const Point& other) const {
            double val = 0;
            for (size_t i = 0; i < coor.size(); ++i) {
                val += coor[i] * other.coor[i];
            }
            return val;
        }

        double norm() const {
            return sqrt(dot(*this));
        }

        void print() const {
            cout << "(";
            for (size_t i = 0; i < coor.size(); ++i) {
                cout << coor[i];
                if (i < coor.size() - 1) cout << ", ";
            }
            cout << ")";
        }
    };

    class QuadraticFunction {
    private:
        vector<vector<double>> A;
        vector<double> b;
        double c;

    public:
        QuadraticFunction(const vector<vector<double>>& A_, const vector<double>& b_, double c_) : A(A_), b(b_), c(c_) {}

        double operator()(const Point& x) const {
            return 0.5 * multiplyAx(x).dot(x) + Point(b).dot(x) + c;
        }
        Point multiplyAx(const Point& x) const {
            vector<double> result(A.size(), 0.0);
            for (size_t i = 0; i < A.size(); ++i) {
                for (size_t j = 0; j < A[i].size(); ++j) {
                    result[i] += A[i][j] * x[j];
                }
            }
            return Point(result);
        }

        Point gradient(const Point& x) const {
            return multiplyAx(x) + Point(b);
        }

        const vector<vector<double>>& getA() const { return A; }
        const vector<double>& getB() const { return b; }
        double getC() const { return c; }
    };


    Point fletcher_reeves_exactly(const QuadraticFunction& f, const Point& x0, double eps = 0.01) {
        Point x = x0;
        Point grad = f.gradient(x);
        Point d = grad * (-1); 

        int k = 0;

        while (grad.norm() > eps) {
            
            Point Ad = f.multiplyAx(d);
            double alpha = grad.dot(grad) / d.dot(Ad);

            x = x + d * alpha;
            Point grad_prev = grad;

            grad = f.gradient(x);

            double omega = grad.dot(grad) / grad_prev.dot(grad_prev);


            d = grad * (-1) + d * omega;

            k++;

        }
		std::cout << "$" << k << "$$";
        return x;
    }

    Point fletcher_reeves_approx(const QuadraticFunction& f, const Point& x0, double eps = 0.01) {
        Point x = x0;
        Point grad = f.gradient(x);
        Point d = grad * (-1); 

        int k = 0;

        while (grad.norm() > eps) {

            double alpha = golden_ratio([&](double alpha) -> double {return f(x + d * alpha);}, 0, 1, eps);

            x = x + d * alpha;

            Point grad_prev = grad;

            grad = f.gradient(x);

            double beta = grad.dot(grad) / grad_prev.dot(grad_prev);

            d = grad * (-1) + d * beta;

            k++;
        }

		std::cout << "$" << k << "$$";
        return x;
    }

    void test() {
        vector<vector<double>> A = {
            {1, 1, 1},
            {1, 1.5, 1},
            {1, 1, 2.5}
        };
        vector<double> b = { 1, -2, -3 };
        double c = 7;
        QuadraticFunction f(A, b, c);
        Point x0(0, 0, 0);
  
        cout << "\nEXACTLy FOR ARE YOU LOOKING FOR" << endl;
        Point result_exact = fletcher_reeves_exactly(f, x0, EPS);
        cout << "ans:";
        result_exact.print();

        cout << "\n ~~~~~~~~~~~~~~~~~~~~~~~~~" << endl;
        Point result_approx = fletcher_reeves_approx(f, x0, EPS);
        cout << "ans: ";
        result_approx.print();
    }
}


int main() {
	//std::cout << TAU;
	//lab_2::test_golden();
	//lab_2::test_fib();
	//lab_3::test_dif_params();
	//std::cout << std::endl;
	//lab_4::test();
	lab_5::test();  

	return 0;
}
