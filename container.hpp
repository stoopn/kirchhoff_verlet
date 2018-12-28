#include <vector>
#include <utility>


template<class C1, class C2>
class wirecontainer {
    std::vector<C1> nodes;
    std::vector<C2> edges;
public:
    
    class wireiterator;
    friend class wireiterator;
    
    class wireiterator : public std::iterator<std::random_access_iterator_tag,C1,std::ptrdiff_t> {
    public:
	typedef  typename std::vector<C1>::iterator C1Iterator;
	typedef  typename std::vector<C2>::iterator C2Iterator;
	typename std::vector<C1>* lst1;
	typename std::vector<C2>* lst2;
	C1Iterator it1;
	C2Iterator it2;
	wireiterator(std::vector<C1>& ls1, std::vector<C2>& ls2, const typename std::vector<C1>::iterator& i1, const typename std::vector<C2>::iterator& i2)
	    : lst1(&ls1), lst2(&ls2), it1(i1), it2(i2) {}

	bool operator==(const wireiterator& x) const {
	    return (it1 == x.it1 && it2==x.it2);
	}
	bool operator!=(const wireiterator& x) const {
	    return (it1 != x.it1 || it2!=x.it2);
	}

	typename std::vector<C1>::reference operator*() const {
	    return *it1;
	}
	
	wireiterator& operator++() {
	    ++it1;
	    ++it2;
	    return *this;
	}
	wireiterator operator+(int n) const {
	    wireiterator tmp = *this;
	    tmp.it1+=n;
	    tmp.it2+=n;
	    return tmp;
	}
	wireiterator operator-(int n) const {
	    wireiterator tmp = *this;
	    tmp.it1-=n;
	    tmp.it2-=n;
	    return tmp;
	  
	}
	std::ptrdiff_t operator-(wireiterator& i1) const {
	   
	    return it1-i1.it1;
	  
	}
	wireiterator& operator-=(int n) {
	    it1-=n;
	    it2-=n;
	    return *this;
	}
	wireiterator& operator+=(int n) {
	    it1+=n;
	    it2+=n;
	    return *this;
	}
	wireiterator operator++(int) {
	    wireiterator tmp = *this;
	    ++*this;
	    return tmp;
	}
	wireiterator& operator--() {
	    --it1;
	    --it2;
	    return *this;
	}
	wireiterator operator--(int) {
	    wireiterator tmp = *this;
	    --*this; 
	    return tmp;
	}
	C1* operator->() {
	    return &(*it1);
	}
	typename std::vector<C2>::reference left_edge() {
	    return *(it2-1);
	}
	typename std::vector<C2>::reference right_edge() {
	    return *it2;
	}
	typename std::vector<C1>::reference left_node() {
	    return *(it1-1);	
	}
	typename std::vector<C1>::reference right_node() {
	    return *(it1+1);
	}
	bool has_right_node() {
	    return (it1+1 < lst1->end());
	}
	bool has_left_node() {
	    return (it1-1 >= lst1->begin());
	}
	bool operator<(const wireiterator& x) const {
	    return (it1 < x.it1);
	}
	bool operator>(const wireiterator& x) const {
	    return (it1 > x.it1);
	}
	typename std::vector<C1> get_nodes() {
	    return nodes;
	}
	typename std::vector<C2> get_edges() {
	    return edges;
	}

    };

    wireiterator begin() {
	return wireiterator(nodes, edges, nodes.begin(), edges.begin());
    }
    wireiterator end() {
	return wireiterator(nodes, edges, nodes.end(), edges.end());
    }
    typename std::vector<C1>::reference operator[](int index) {
	return nodes[index];
    }
    typename std::vector<C2>::reference get_edge(int index) {
	return edges[index];
    }
    typename std::pair<C1,C2> get_pair(int index) {
	return std::pair<C1,C2>(nodes[index],edges[index]);
    }
    wireiterator insert(wireiterator pos, C1& n, C2& e) {
	typename std::vector<C1>::iterator i1 = nodes.insert(pos.it1,n);
	typename std::vector<C2>::iterator i2 = edges.insert(pos.it2,e);
	return wireiterator(nodes, edges, i1, i2);
    }
    wireiterator erase(wireiterator pos) {
	typename std::vector<C1>::iterator i1 = nodes.erase(pos.it1);
	typename std::vector<C2>::iterator i2 = edges.erase(pos.it2);
	return wireiterator(nodes, edges, i1, i2);
    }
    void resize(int n, C1 t1=C1(),C2 t2=C2()) {
	nodes.resize(n,t1);
	edges.resize(n,t2);
    }
    
    int size() { return nodes.size(); }
    bool empty() { return nodes.empty(); }
    typename std::vector<C1> get_nodes() const {return nodes;}
    typename std::vector<C2> get_edges() const {return edges;}
};
