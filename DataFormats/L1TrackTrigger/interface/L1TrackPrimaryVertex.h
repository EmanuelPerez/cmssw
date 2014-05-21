#ifndef L1TkTrigger_L1TrackPrimaryVertex_h
#define L1TkTrigger_L1TrackPrimaryVertex_h


// Nov 12, 2013
// First version of a class for L1-zvertex

class L1TrackPrimaryVertex {

 public:

 L1TrackPrimaryVertex() : zvertex(-999), sum(-999) {}

 ~L1TrackPrimaryVertex() { }

 L1TrackPrimaryVertex(float z, float s) : zvertex(z), sum(s) { }


    float getZvertex() const { return zvertex ; } 
    float getSum() const { return sum ; }

 private:
   float zvertex;
   float sum;

};

#include <vector>

typedef std::vector<L1TrackPrimaryVertex> L1TrackPrimaryVertexCollection ;


#endif
