#include "fusedanova.h"

// Constructor/Destructor
node:: node() {};
node::~node() {};
node::node(int n        ,
           int label_   , 
           double beta_ ,
           double slope_, 
           int weight_  ) 
    : lambda(0.0)    ,
      beta(beta_)    ,
      slope(slope_)  ,
      range(n)       ,
      weight(weight_),
      size(1)        , 
      label(label_)  ,
      idown(label_)  ,
      isplit(label_) ,
      iup(label_)    ,
      parent1(label_),
      parent2(label_),
      active(true) 
  {
    if (label == 0  ) down = -1; else  down = label - 1;
    if (label == n-1) up   = -1; else  up   = label + 1;
  } 
;

// + operator to fuse two nodes
node node::operator+ (const node& node_) {

  int weight_ = this->weight + node_.weight ;
  double lambda_ = (this->beta - node_.beta - this->slope * this->lambda + node_.slope * node_.lambda) / (node_.slope - this->slope);
  double slope_  = (this->weight * this->slope + node_.weight * node_.slope) / weight_ ;
  double beta_   = this->beta + (lambda_ - this->lambda) * this->slope ;
  node result(this->range, 0, beta_, slope_, weight_) ;
  
  result.lambda = lambda_ ;
  result.size  = this->size + node_.size;
  result.parent1 = this->label;
  result.parent2 = node_.label;
  if(this->idown < node_.idown) {
    result.isplit = this->iup   ;
    result.idown  = this->idown ;
    result.iup    = node_.iup   ;
    result.down   = this->down  ;
    result.up     = node_.up    ;
  } else {
    result.isplit = node_.iup   ; 
    result.idown  = node_.idown ;
    result.iup    = this->iup   ;
    result.down   = node_.down  ;
    result.up     = this->up    ;
  }
  
  return(result);  
};
  
Fusion::Fusion(node *node1_, node *node2_) 
  : node1(node1_),  node2(node2_) 
  {
    lambda = (node1->beta - node2->beta - node1->slope * node1->lambda + node2->slope * node2->lambda) / (node2->slope - node1->slope) ;
  } 
;

