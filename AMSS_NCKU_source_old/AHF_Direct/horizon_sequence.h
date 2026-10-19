#ifndef HORIZON_SEQUENCE_H
#define HORIZON_SEQUENCE_H
namespace AHFinderDirect
{
	class horizon_sequence
	{
	public:
		int N_horizons() const { return N_horizons_; }

		int my_N_horizons() const { return my_N_horizons_; }

		bool has_genuine_horizons() const { return my_N_horizons_ > 0; }

		bool is_dummy() const { return posn_is_dummy(posn_); }
		bool is_genuine() const { return posn_is_genuine(posn_); }

		bool is_next_genuine() const
		{
			return posn_is_genuine(next_posn(posn_));
		}

		int dummy_number() const { return is_genuine() ? 0 : -posn_; }

		int get_hn() const
		{
			return posn_is_genuine(posn_) ? my_hn_[posn_] : 0;
		}

		bool is_hn_genuine(int hn) const;

		int init_hn()
		{
			posn_ = (my_N_horizons_ == 0) ? -1 : 0;
			return get_hn();
		}

		int next_hn()
		{
			posn_ = next_posn(posn_);
			return get_hn();
		}

		horizon_sequence(int N_horizons);
		~horizon_sequence();

		int append_hn(int hn);

	private:
		bool posn_is_genuine(int pos) const
		{
			return (pos >= 0) && (pos < my_N_horizons_);
		}
		bool posn_is_dummy(int pos) const
		{
			return !posn_is_genuine(pos);
		}

		int next_posn(int pos) const;

	private:
		const int N_horizons_;
		int my_N_horizons_;

		int posn_;

		int *my_hn_;
	};

	//******************************************************************************

} // namespace AHFinderDirect
#endif /* HORIZON_SEQUENCE_H */
