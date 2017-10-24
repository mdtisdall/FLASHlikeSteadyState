#ifndef SequenceEventVectors_h
#define SequenceEventVectors_h

template <typename vec_type_>
class SequenceEventVectors {
  public:
    typedef vec_type_ vec_type;
    typedef typename vec_type::value_type value_type;

    SequenceEventVectors(const value_type t1_ms) :
      invT1_per_ms(value_type(1.0) / t1_ms)
    {}

    class Excitation {
      public:
        Excitation(
          const value_type flipAngle_deg,
          const value_type postExcitationDelay_ms) :
          cosFlipAngle(std::cos(flipAngle_deg * M_PI / value_type(180.0)))
         ,postExcitationDelay_ms(postExcitationDelay_ms)
        {}

        virtual void operator()(SequenceEventVectors<vec_type> *vecs) const {
          value_type scale =
            std::exp(-postExcitationDelay_ms * vecs->invT1_per_ms);
          vecs->a.push_back(-scale * cosFlipAngle);
          vecs->b.push_back(value_type(1.0) - scale);
        }
        
        virtual void operator()(
          typename vec_type::iterator *longIt,
          vec_type*) const
        {
          (*longIt)++;
        }

      protected:
        const value_type cosFlipAngle;
        const value_type postExcitationDelay_ms;
    };

    class ExcitationWithReadout : public Excitation
    {
      public:
        ExcitationWithReadout(
          const value_type flipAngle_deg,
          const value_type postExcitationDelay_ms) :
          Excitation(flipAngle_deg, postExcitationDelay_ms)
         ,sinFlipAngle(std::sin(flipAngle_deg * M_PI / value_type(180.0)))
        {}
       
        using Excitation::operator();

        virtual void operator()(
          typename vec_type::iterator *longIt,
          vec_type *transverse) const
        {
          transverse->push_back(**longIt * sinFlipAngle);
          (*longIt)++;
        }

      protected:
        const value_type sinFlipAngle;
    };
    
    class Relaxation {
      public:
        Relaxation(
          const value_type duration_ms) :
          duration_ms(duration_ms)
        {}

        void operator()(SequenceEventVectors<vec_type> *vecs) const {
          value_type scale =
            std::exp(-duration_ms * vecs->invT1_per_ms);
          vecs->a.push_back(-scale);
          vecs->b.push_back(value_type(1.0) - scale);
        }
        
        void operator()(
          typename vec_type::iterator *longIt,
          vec_type*) const {
          (*longIt)++;
        }

      protected:
        const value_type duration_ms;
    };

    const value_type invT1_per_ms;
    vec_type a;
    vec_type b;

};

#endif
