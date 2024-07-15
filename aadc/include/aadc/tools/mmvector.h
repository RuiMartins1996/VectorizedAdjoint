#ifndef AADC_MMVECTOR_H
#define AADC_MMVECTOR_H

#include <vector>

#ifdef AADC_USE_BOOST_ALIGNED_ALLOCATOR
#include <boost/align/aligned_allocator.hpp>
namespace aadc {

template <class StoreType, class mmType = StoreType>
using mmVector = std::vector<StoreType, boost::alignment::aligned_allocator<mmType, sizeof(mmType)>>; // mmVector<mmType>;

};

#else // AADC_USE_BOOST_ALIGNED_ALLOCATOR

namespace aadc {

template <typename T, std::size_t alignment>
class aligned_allocator
{
public:
	typedef T * pointer;
	typedef const T * const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	typedef T value_type;
	typedef std::size_t size_type;
	typedef std::ptrdiff_t difference_type;
public:
	template <typename U>
	struct rebind
	{
		typedef aligned_allocator<U, alignment> other;
	};

public:
	aligned_allocator() = default;
	aligned_allocator(const aligned_allocator&) = default;

	template <typename U> 
	aligned_allocator(const aligned_allocator<U, alignment>&) {};
	~aligned_allocator() = default;

	T * address(T& r) const
	{
		return &r;
	}

	const T * address(const T& s) const
	{
		return &s;
	}

	std::size_t max_size() const
	{
		return (static_cast<std::size_t>(0) - static_cast<std::size_t>(1)) / sizeof(T);
	}

	bool operator!=(const aligned_allocator& other) const
	{
		return !(*this == other);
	}

	void construct(T * const p, const T& t) const
	{
		void * const pv = static_cast<void *>(p);

		new (pv) T(t);
	}

	void destroy(T * const p) const
	{
		p->~T();
	}

	bool operator==(const aligned_allocator& other) const
	{
		return true;
	}

	T * allocate(const std::size_t n) const
	{
		if (n == 0) return NULL;
		if (n > max_size()) throw std::bad_alloc();

		void * const pv = _mm_malloc(n * sizeof(T), alignment);
		if (pv == NULL) throw std::bad_alloc();

		return static_cast<T *>(pv);
	}

	void deallocate(T * const p, const std::size_t n) const
	{
		_mm_free(p);
	}

	template <typename U>
	T * allocate(const std::size_t n, const U *) const
	{
		return allocate(n);
	}

private:
	aligned_allocator& operator=(const aligned_allocator&);
};

template <class StoreType, class mmType = StoreType>
using mmVector = std::vector<StoreType, aligned_allocator<mmType, sizeof(mmType)>>;

};

#endif // AADC_USE_BOOST_ALIGNED_ALLOCATOR
#endif //AADC_MMVECTOR_H
