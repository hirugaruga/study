#include <iostream>

template<typename T>
class Stack {
 public:
	Stack();
	void push(T val);
	T pop();
	T back() const;
	int size() const;
	void clear();
	bool empty();
	T min() const;
	T max() const;
 private:
	struct Node {
		T val;
		Node* next;
		T min;
		T max;
	};
	Node* top_;
	int size_;
};

int main() {
	Stack<int> s;
	s.push(1);
	s.push(8);
	s.push(7);
	s.push(6);
	s.push(5);
	s.push(4);
	std::cout << s.min() << std::endl;
	s.pop();
	std::cout << s.min() << std::endl;
	s.pop();
	std::cout << s.min() << std::endl;
	s.pop();
	std::cout << s.min() << std::endl;
	s.pop();
	std::cout << s.min() << std::endl;
	s.pop();

}
/****************************Stack*********************************************/
template<typename T>
Stack<T>::Stack() {
	top_ = nullptr;
	size_ = 0;
}

template<typename T>
void Stack<T>::push(T val) {
	Node* tmp = new Node;
	tmp->val = val;
	tmp->next = top_;
	top_ = tmp;
	++size_;
	if (size_ <= 1) {
		top_->min = val;
		top_->max = val;
	} else {
		top_->min = (top_->next->min > val) ? val : top_->next->min;
		top_->max = (top_->next->max < val) ? val : top_->next->max;
	}
}

template<typename T>
T Stack<T>::pop() {
	if (size_ == 0) {
		std::cerr << "Stack::pop : empty stack \n";
	} else {
		Node* tmp = top_;
		top_ = top_->next;
		--size_;
		int val = tmp->val;
		return val;
		delete tmp;
	}
}

template<typename T>
T Stack<T>::back() const {
	return top_->val;
}

template<typename T>
int Stack<T>::size() const {
	return size_;
}

template<typename T>
void Stack<T>::clear() {
	while (size_ != 0) {
		Node* tmp = top_;
		top_ = top_->next;
		--size_;
		delete tmp;
	}
	top_ = nullptr;
}

template<typename T>
bool Stack<T>::empty() {
	if (size_ == 0) {
		return true;
	} else {
		return false;
	}
}

template<typename T>
T Stack<T>::min() const {
	return top_->min;
}

template<typename T>
T Stack<T>::max() const {
	return top_->max;
}