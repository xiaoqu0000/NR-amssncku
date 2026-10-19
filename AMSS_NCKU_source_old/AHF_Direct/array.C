#include <assert.h>
#include <stddef.h> // NULL
#include <stdlib.h> // size_t

#include "cctk.h"

#include "stdc.h"
#include "util.h"
#include "array.h"

namespace AHFinderDirect
{
	namespace jtutil
	{

		template <typename T>
		array1d<T>::array1d(int min_i_in, int max_i_in,
							T *array_in /* = NULL */,
							int stride_i_in /* = 0 */)
			: array_(array_in),
			  offset_(0), // temp value, changed below
			  stride_i_(stride_i_in),
			  min_i_(min_i_in), max_i_(max_i_in),
			  we_own_array_(array_in == NULL)
		{
			if (stride_i_ == 0)
				then stride_i_ = 1;

			// must use unchecked subscripting here since setup isn't done yet
			offset_ = -subscript_unchecked(min_i_); // RHS uses offset_ = 0
			assert(subscript_unchecked(min_i_) == 0);
			max_subscript_ = subscript_unchecked(max_i_);

			if (we_own_array_)
				then
				{
					// allocate it
					const int N_allocate = N_i();
					array_ = new T[N_allocate];
				}

			// explicitly initialize array (new[] *doesn't* do this automagically)
			for (int i = min_i(); i <= max_i(); ++i)
			{
				operator()(i) = T(0);
			}
		}

		//
		// This function destroys an  array1d  object.
		//
		template <typename T>
		array1d<T>::~array1d()
		{
			if (we_own_array_)
				then delete[] array_;
		}

		//
		// This function constructs an  array2d  object.
		//
		template <typename T>
		array2d<T>::array2d(int min_i_in, int max_i_in,
							int min_j_in, int max_j_in,
							T *array_in /* = NULL */,
							int stride_i_in /* = 0 */, int stride_j_in /* = 0 */)
			: array_(array_in),
			  offset_(0), // temp value, changed below
			  stride_i_(stride_i_in), stride_j_(stride_j_in),
			  min_i_(min_i_in), max_i_(max_i_in),
			  min_j_(min_j_in), max_j_(max_j_in),
			  we_own_array_(array_in == NULL)
		{
			if (stride_j_ == 0)
				then stride_j_ = 1;
			if (stride_i_ == 0)
				then stride_i_ = N_j();

			// must use unchecked subscripting here since setup isn't done yet
			offset_ = -subscript_unchecked(min_i_, min_j_); // RHS uses offset_ = 0
			assert(subscript_unchecked(min_i_, min_j_) == 0);
			max_subscript_ = subscript_unchecked(max_i_, max_j_);

			if (we_own_array_)
				then
				{
					// allocate it
					const int N_allocate = N_i() * N_j();
					array_ = new T[N_allocate];
				}

			// explicitly initialize array (new[] *doesn't* do this automagically)
			for (int i = min_i(); i <= max_i(); ++i)
			{
				for (int j = min_j(); j <= max_j(); ++j)
				{
					operator()(i, j) = T(0);
				}
			}
		}

		//
		// This function destroys an  array2d  object.
		//
		template <typename T>
		array2d<T>::~array2d()
		{
			if (we_own_array_)
				then delete[] array_;
		}

		//
		// This function constructs an  array3d  object.
		//
		template <typename T>
		array3d<T>::array3d(int min_i_in, int max_i_in,
							int min_j_in, int max_j_in,
							int min_k_in, int max_k_in,
							T *array_in /* = NULL */,
							int stride_i_in /* = 0 */, int stride_j_in /* = 0 */,
							int stride_k_in /* = 0 */)
			: array_(array_in),
			  offset_(0), // temp value, changed below
			  stride_i_(stride_i_in), stride_j_(stride_j_in),
			  stride_k_(stride_k_in),
			  min_i_(min_i_in), max_i_(max_i_in),
			  min_j_(min_j_in), max_j_(max_j_in),
			  min_k_(min_k_in), max_k_(max_k_in),
			  we_own_array_(array_in == NULL)
		{
			if (stride_k_ == 0)
				then stride_k_ = 1;
			if (stride_j_ == 0)
				then stride_j_ = N_k();
			if (stride_i_ == 0)
				then stride_i_ = N_j() * N_k();

			// must use unchecked subscripting here since setup isn't done yet
			offset_ = -subscript_unchecked(min_i_, min_j_, min_k_); // RHS uses offset_ = 0
			assert(subscript_unchecked(min_i_, min_j_, min_k_) == 0);
			max_subscript_ = subscript_unchecked(max_i_, max_j_, max_k_);

			if (we_own_array_)
				then
				{
					// allocate it
					const int N_allocate = N_i() * N_j() * N_k();
					array_ = new T[N_allocate];
				}

			// explicitly initialize array (new[] *doesn't* do this automagically)
			for (int i = min_i(); i <= max_i(); ++i)
			{
				for (int j = min_j(); j <= max_j(); ++j)
				{
					for (int k = min_k(); k <= max_k(); ++k)
					{
						operator()(i, j, k) = T(0);
					}
				}
			}
		}
		//
		// This function destroys an  array3d  object.
		//
		template <typename T>
		array3d<T>::~array3d()
		{
			if (we_own_array_)
				then delete[] array_;
		}

		template class array1d<int>;

		// FIXME: we shouldn't have to instantiate these both, the const one
		//	  is actually trivially derivable from the non-const one. :(
		template class array1d<void *>;
		template class array1d<const void *>;

		template class array1d<CCTK_REAL>;
		template class array2d<CCTK_INT>;
		template class array2d<CCTK_REAL>;
		template class array3d<CCTK_REAL>;

	} // namespace jtutil
} // namespace AHFinderDirect
