#ifndef MPRAGE_h
#define MPRAGE_h

template <typename SequenceEventVectors>
class MPRAGE {
  public:
    typedef typename SequenceEventVectors::vec_type vec_type;
    typedef typename SequenceEventVectors::value_type value_type;

    MPRAGE(
      const value_type flipAngle_deg,
      const value_type td1_ms,
      const value_type td2_ms,
      const unsigned int innerSteps,
      const value_type flashTR_ms
    ) :
      innerSteps(innerSteps)
     ,excitationWithReadout(flipAngle_deg, flashTR_ms)
     ,inversion(180, td1_ms)
     ,td2Relaxation(td2_ms)
    {}

    void operator()(SequenceEventVectors *vecs)
    {
      inversion(vecs);

      for(unsigned int innerStep = 0; innerStep < innerSteps; innerStep++) {
         excitationWithReadout(vecs); 
      }
      
      td2Relaxation(vecs);
    }
  protected:
    const unsigned int innerSteps;
    typename SequenceEventVectors::ExcitationWithReadout excitationWithReadout;
    typename SequenceEventVectors::Excitation inversion;
    typename SequenceEventVectors::Relaxation td2Relaxation;
};

#endif
