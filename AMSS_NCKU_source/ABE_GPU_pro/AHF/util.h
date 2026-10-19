#ifndef AHFINDERDIRECT__UTIL_HH
#define AHFINDERDIRECT__UTIL_HH
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <strstream>
#include <cmath>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <string.h>
#include <math.h>
#endif

#define PI M_PI

namespace AHFinderDirect
{
	namespace jtutil
	{
		inline int how_many_in_range(int low, int high) { return high - low + 1; }

		inline int is_even(int i) { return !(i & 0x1); }
		inline int is_odd(int i) { return (i & 0x1); }

		template <typename T>
		inline T tmin(T x, T y) { return (x < y) ? x : y; }
		template <typename T>
		inline T tmax(T x, T y) { return (x > y) ? x : y; }
		template <typename T>
		inline T abs(T x) { return (x > 0) ? x : -x; }

		template <typename T>
		inline T pow2(T x) { return x * x; }
		template <typename T>
		inline T pow3(T x) { return x * x * x; }
		template <typename T>
		inline T pow4(T x) { return pow2(pow2(x)); }

		template <typename fp_t>
		inline fp_t degrees_of_radians(fp_t radians) { return (180.0 / PI) * radians; }
		template <typename fp_t>
		inline fp_t radians_of_degrees(fp_t degrees) { return (PI / 180.0) * degrees; }

		// in miscfp.cc
		//-----------------------------------------------------
		double signum(double x);
		double hypot3(double x, double y, double z);
		double arctan_xy(double x, double y);

		double modulo_reduce(double x, double xmod, double xmin, double xmax);

		template <typename fp_t>
		void zero_C_array(int N, fp_t array[]);

		// in error_exit.cc
		// ------------------------------------------------------
		int error_exit(int msg_level, const char *format, ...);

		// in norm.cc
		//
		template <typename fp_t>
		class norm
		{
		public:
			// get norms etc
			fp_t mean() const;
			fp_t two_norm() const; // sqrt(sum x_i^2)
			fp_t rms_norm() const; // sqrt(average of x_i^2)
			fp_t infinity_norm() const { return max_abs_value_; }

			fp_t max_abs_value() const { return max_abs_value_; }
			fp_t min_abs_value() const { return min_abs_value_; }

			fp_t max_value() const { return max_value_; }
			fp_t min_value() const { return min_value_; }

			// specify data point
			void data(fp_t x);

			// have any data points been specified?
			bool is_empty() const { return N_ == 0; }
			bool is_nonempty() const { return N_ > 0; }

			// reset ==> just like newly-constructed object
			void reset();

			// constructor, destructor
			// ... compiler-generated no-op destructor is ok
			norm();

		private:
			// we forbid copying and passing by value
			// by declaring the copy constructor and assignment operator
			// private, but never defining them
			norm(const norm &rhs);
			norm &operator=(const norm &rhs);

		private:
			long N_;			 // # of data points
			fp_t sum_;			 // sum(data)
			fp_t sum2_;			 // sum(data^2)
			fp_t max_abs_value_; // max |data|
			fp_t min_abs_value_; // min |data|
			fp_t max_value_;	 // max data
			fp_t min_value_;	 // min data
		};

		// in fuzzy.cc
		template <typename fp_t>
		class fuzzy
		{
		public:
			// comparison tolerance (may be modified by user code if needed)
			static fp_t get_tolerance() { return tolerance_; }
			static void set_tolerance(fp_t new_tolerance)
			{
				tolerance_ = new_tolerance;
			}

			// fuzzy commparisons
			static bool EQ(fp_t x, fp_t y);
			static bool NE(fp_t x, fp_t y) { return !EQ(x, y); }
			static bool LT(fp_t x, fp_t y) { return EQ(x, y) ? false : (x < y); }
			static bool LE(fp_t x, fp_t y) { return EQ(x, y) ? true : (x < y); }
			static bool GT(fp_t x, fp_t y) { return EQ(x, y) ? false : (x > y); }
			static bool GE(fp_t x, fp_t y) { return EQ(x, y) ? true : (x > y); }

			static bool is_integer(fp_t x); // is x fuzzily an integer?
			static int floor(fp_t x);		// round x fuzzily down to integer
			static int ceiling(fp_t x);		// round x fuzzily up to integer

		private:
			// comparison tolerance
			// ... must be explicitly initialized when instantiating
			//     for a new <fp_t> type, see "fuzzy.cc" for details/examples
			static fp_t tolerance_;
		};

		// in round.cc
		template <typename fp_t>
		class round
		{
		public:
			static int to_integer(fp_t x); // round to nearest integer

			static int floor(fp_t x);	// round down to integer
			static int ceiling(fp_t x); // round up to integer
		};

	} // namespace jtutil
} // namespace AHFinderDirect

#endif /* AHFINDERDIRECT__UTIL_HH */
