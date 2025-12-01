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

#include <iostream>
#include <optional>
#include <random>
#include <cmath>

#include <stack>


namespace problem_one {
	//задание 1.1. Точный алгоритм для разности
	double smth(int a, int b) {
		if (a == b && a % 2 == 0) return 0;
		if (a == b && a % 1 == 0) return 1;
		return 0;

	}
	//задание 1.2. Приблеженный алгоритм для разности
	// Методом кавальери или рандомом ;/

	double diff_black_white(double a, double b) {
		int n = 1e6;
		std::random_device rd;
		std::mt19937_64 gen(rd());
		std::uniform_real_distribution<double> dist_x(0, a);
		std::uniform_real_distribution<double> dist_y(0, b);
		long long black = 0;
		long long white = 0;

		for (int i = 0; i < n; i++) {
			double x = dist_x(gen);
			double y = dist_y(gen);
			if (y <= b * (1 - x / a)) {
				int rx = (int)x;
				int ry = (int)y;
				if ((rx + ry) % 2 == 0) ++black;
				else ++white;
			}
			else n--;
		}
		double s = a * b / 2;
		double black_s = s * black / n * 2;
		double white_s = s * white / n * 2;
		return black_s - white_s;
		}



	}
//задание 2.1. Дерево поиск  

namespace problem_two {

	template <typename TKey, typename TData>
	class binary_search_tree
	{
		struct node {
			TKey key; TData data; node* left, right; size_t size;
			node(const TKey& k, const TData& d) : key(k), data(d), left(nullptr), right(nullptr) {};
		}
		*root = nullptr;
		size_t count = 0;

		size_t get_size(node* n) const { return n ? n->size : 0; }


	public:


		binary_search_tree() {};
		
		~binary_search_tree();

		std::optional<TData> find(const TKey& key);

		void insert(const TKey& key, const TData& data);

		std::optional<TKey> findnext(const TKey& key);

		void remove(const TKey& key);

		//void remove_alt(const TKey& key);

	};

	template<typename TKey, typename TData>
	std::optional<TData> binary_search_tree<TKey, TData>::find(const TKey& key)
	{
		node* current = root;
		while (current != nullptr && current->key != key) {
			current = current->key > key ? current->left: current->right;
		}
		return (current != nullptr && current->key == key) ? current->data : std::nullopt;
	}

	template<typename TKey, typename TData>
	void binary_search_tree<TKey, TData>::insert(const TKey& key, const TData& data)
	{
		node* current = root;
		node* parent = nullptr;
		while (current != nullptr && current->key != key) {
			parent = current;
			current = current->key > key ? current->left : current->right;
		}
		if (current != nullptr && current->key == key) current->data = data;
		else {
			node* tmp = new node{ key, data };
			if (parent == nullptr) root = tmp;
			else if (key < parent->key) parent->left = tmp;
			else parent->right = tmp;
			count++;
		}
	}

	template<typename TKey, typename TData>
	std::optional<TKey> binary_search_tree<TKey, TData>::findnext(const TKey& key)
	{
		node* current = root;
		node* last_left = nullptr;
		while (current != nullptr && current->key != key) {
			if (current->key > key) {
				last_left = current;
				current = current->left;
			}
			else current = current->right;
		}
		if (current->right != nullptr) {
			current = current->right;
			while (current->left != nullptr) current = current->left;
			return current->key;
		}
		else return last_left != nullptr ? last_left->key: std::nullopt;
	}

	template<typename TKey, typename TData> 
	void binary_search_tree<TKey, TData>::remove(const TKey& key)
	{
		node* current = root;
		node* parent = nullptr;
		while (current != nullptr && current->key != key) {
			parent = current;
			current = current->key > key ? current->left : current->right;
		}
		if (current == nullptr) return;
		if (current->left == nullptr && current->right == nullptr) {
			if (parent == nullptr) root = nullptr;
			else if (parent->left == current) parent->left = nullptr;
			else parent->right = nullptr;
		} else if (current->left == nullptr || current->right == nullptr) {
			node* child = current->left != nullptr ? current->left : current->right;
			if (parent == nullptr) root = child;
			else if (parent->left == current) parent->left = child;
			else parent->right = child;
		} else {
			node* A = current->right;
			node* A_parent = current;
			while (A->left != nullptr) {
				A_parent = A;
				A = A->left;
			}
			current->key = A->key;
			current->data = A->data;

			if (A_parent->left == A) A_parent->left = A->right;
			else A_parent->right = A->right;
			delete A;
			count--;
			return;
		}
		delete current;
		count--;

	}

	template<typename TKey, typename TData>
	binary_search_tree<TKey, TData>::~binary_search_tree()
	{
		std::stack<node*> S;
		node* current = root;
		while (current != nullptr || !S.empty()) {
			if (current != nullptr) {
				S.push(current);
				current = current->left;
			} else {
				current = S.top()->right;
				delete S.top();
				S.pop();
			}
		}
	}
}


int main() {
	std::cout << problem_one::diff_black_white(3, 3) << std::endl;

}


