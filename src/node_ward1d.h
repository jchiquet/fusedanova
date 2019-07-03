#ifndef _node_ward1d_H
#define _node_ward1d_H

#include <cmath>

class node_ward1d {
public:
  double lambda ; 
  double xplus  ;
  double x2plus ;
  int K         ;
  int weight    ;
  int size      ;
  int label     ;
  int idown     ;
  int isplit    ;
  int iup       ;
  int down      ;
  int up        ;
  int parent1   ;
  int parent2   ;
  bool active   ;

  // Constructors/Destructor
   node_ward1d() ;
  ~node_ward1d() ;
  node_ward1d(int n, int label_, double xplus_, double x2plus_, int weight_) ;

  // Basic methods for acces
  bool has_down () const {return (down != -1);} ;
  bool has_up   () const {return (up   != -1);} ;

  // Define overloaded + for fusing two nodes
  node_ward1d operator+ (const node_ward1d& node_) ;
};

class Fusion_ward1d {

public:
  node_ward1d *node1 ;
  node_ward1d *node2 ;
  double lambda ;

  // Constructor
  Fusion_ward1d(node_ward1d *node1_, node_ward1d *node2_) ; 
  // Getter
  int label1() {return node1->label ;}  
  int label2() {return node2->label ;}  
  double get_lambda() const {return lambda ;}
  bool is_active() const {return(node1->active & node2->active);}
};

class UpcomingFusions_ward1d {
public:
  double operator() (const Fusion_ward1d& r1, const Fusion_ward1d& r2) {
    return r1.get_lambda() >= r2.get_lambda();
  }
};

#endif        
