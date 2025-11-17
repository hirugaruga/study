#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <set>

std::string MARKERS = "+-*/^()";

namespace problem_one {
	enum class lexem_t {
		VAR, CONST, OPERATOR, FUNC, OP_BR, CL_BR
	};


	struct lexem
	{
		std::string value;
		lexem_t type;
		
	};

	std::ostream& operator<<(std::ostream& os, lexem_t type) {
		switch (type) {
		case lexem_t::CONST: return os << "CONST"; break;
		case lexem_t::VAR: return os << "VAR"; break;
		case lexem_t::OPERATOR: return os << "OPERATOR"; break;
		case lexem_t::FUNC: return os << "FUNC"; break;
		case lexem_t::OP_BR: return os << "OP_BR"; break;
		case lexem_t::CL_BR: return os << "CL_BR"; break;
		}
	}

	std::ostream& operator<<(std::ostream& os, lexem& l) {
		return os << "{ \"" << l.value << "\" , " << l.type << " }";
	}


	lexem_t type_of_separator(char ch) {
		switch (ch) {
			case '(': return lexem_t::OP_BR; break;
			case ')': return  lexem_t::CL_BR; break;
			default: return lexem_t::OPERATOR;
		}
	}

	std::vector<lexem> token_search(std::string& s) {
		std::vector<lexem> ans;
		//std::vector<std::string> ans;
		auto it = s.begin();
		
		while (it != s.end()) {

			auto end = std::find_first_of(it, s.end(), MARKERS.begin(), MARKERS.end());
			
			if (static_cast<int>(*it) >= 48 && static_cast<int>(*it) <= 57) 
			{
				ans.push_back({ std::string(it, end), lexem_t::CONST });
			}
			else if (*end == '(' && it != end) 
			{
				ans.push_back({ std::string(it, end), lexem_t::FUNC });
			}
			else if (it != end) 
			{
				ans.push_back({ std::string(it, end), lexem_t::VAR });
			}
			ans.push_back({ std::string(end, end + 1), type_of_separator(*end) });

			it = end + 1;
		}
		return ans;
	}

	void test() {
		std::string s = "3.15*abc*foo(x^2/(y-z(x)))-exp(-x)";
		std::vector<lexem> v = token_search(s);
		for (auto elem : v) {
			std::cout << elem << "\n";
		}

	}
}

namespace problem_two {
	struct Point {
		std::vector<double> v;

		Point(std::vector<double> s) : v(s) {};
		Point(double x, double y, double z) {
			v.push_back(x);
			v.push_back(y);
			v.push_back(z);
		};

		double dot() {
			double val = 0;
			for (auto x : v) val += x * x;
			return val;
		}

		size_t size() {
			return v.size();
		}

		double dot(Point& other) {
			double val = 0;
			for (size_t i = 0; i < v.size(); i++) val += v[i] * other[i];
			return val;
		}

		double operator[](size_t ind) {
			return v[ind];
		}


		friend Point operator-(Point& a, Point& b) {
			std::vector<double> res;
			for (size_t i = 0; i < a.size(); i++) res.push_back(a[i] - b[i]);
			return { res };
		}

		//friend Point operator-(Point& a, Point& b)
	};


	long long fac(int n) {
		long long val = 1;
		for (unsigned int i = 1; i < n + 1; i++) val *= i;
		return val;

	}

	long long number_of_uniq_planes(std::vector<Point>& vp) {
		size_t n = vp.size();
		//long long answer = fac(n) / (fac(n - 3) * 6);

		std::set<std::set<size_t>> plates;

		for (size_t i = 0; i < n; ++i) {
			for (size_t j = i + 1; j < n; j++) {
				for (size_t k = j + 1; k < n; ++k) {

					Point x = vp[j] - vp[i];
					Point y = vp[k] - vp[i];
					double m11 = x.dot(x);
					double m12 = x.dot(y);
					double m22 = y.dot(y);
					double det = m11 * m22 - m12 * m12;

					if (std::abs(det) < 1e-9) continue;
					std::set<size_t> temp{ i, j ,k };
					for (size_t m = 0; m < n; m++) {

						if (m == i || m == j || m == k) {
							continue;
						}
						Point z = vp[m] - vp[i];


						double xz = x.dot(z);
						double yz = y.dot(z);




						double q = (m22 * xz - m12 * yz) / det;
						double w = (m11 * yz - m12 * xz) / det;

						bool flag = true;
						for (size_t l = 0; l < x.size(); l++) {
							if (std::abs(z[l] - (q * x[l] + w * y[l])) > 1e-9) { flag = false; break; }
						}
						if (flag) temp.insert(m);
						//if (flag) answer--a;
					}
					if (plates.find(temp) == plates.end()) {
						for (auto elem : temp) {
							std::cout << elem + 1 << ' ';
						}
						std::cout << "\n";
						plates.insert(temp);
					}
					/*for (auto elem : temp) {
						std::cout << elem << ' ';
					}
					std::cout << "\n";
					plates.insert(temp);*/


				}

			}
		}
		return plates.size();
	}


	void test() {
		std::vector<Point> t = {
			{0, 0, 0},
			{0, 0, 1},
			{0, 1, 0},
			{0, 1, 1},
			{ 1, 0, 0 },
			{ 1, 0, 1 },
			{ 1, 1, 0 },
			{ 1, 1, 1 },

		};


		std::cout << number_of_uniq_planes(t);
	
	}

}


int main() {
	//problem_one::test();
	problem_two::test();

}


