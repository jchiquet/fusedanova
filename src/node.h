#ifndef _node_H
#define _node_H

class node {
public:
  double lambda ; 
  double beta   ;
  double slope  ;
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
   node() ;
  ~node() ;
  node(int n, int label_, double beta_, double slope_, int weight_) ;

  // Basic methods for acces
  bool has_down () const {return (down != -1);} ;
  bool has_up   () const {return (up   != -1);} ;

  // Define overloaded + for fusing two nodes
  node operator+ (const node& node_) ;
};

class Fusion {

public:
  node *node1 ;
  node *node2 ;
  double lambda ;

    // Constructor
  Fusion(node *node1_, node *node2_) ; 
  // Getter
  int label1() {return node1->label ;}  
  int label2() {return node2->label ;}  
  double get_lambda() const {return lambda ;}
  bool is_active() const {return(node1->active & node2->active);}
};

class UpcomingFusions {
public:
  double operator() (const Fusion& r1, const Fusion& r2) {
    return r1.get_lambda() >= r2.get_lambda();
  }
};

#endif        
