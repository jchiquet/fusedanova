#include "node_ward1d.h"

// Constructor/Destructor
node_ward1d:: node_ward1d() {};
node_ward1d::~node_ward1d() {};
node_ward1d::node_ward1d(int K_,
           int label_   , 
           double xplus_ ,
           double x2plus_,
           int weight_  ) 
    : lambda(0.0)    ,
      xplus(xplus_)  ,
      x2plus(x2plus_),
      inertia(x2plus_ - pow(xplus_, 2)/weight_),
      K(K_)          ,
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
    if (label == K-1) up   = -1; else  up   = label + 1;
  } 
;

// + operator to fuse two nodes
node_ward1d node_ward1d::operator+ (const node_ward1d& node_) {

  int weight_    = this->weight + node_.weight ;
  double xplus_  = this->xplus  + node_.xplus ;
  double x2plus_ = this->x2plus + node_.x2plus ;

  node_ward1d result(this->K, 0, xplus_, x2plus_, weight_) ;

  result.lambda  =  sqrt(2 *std::max(0.0, result.inertia - this->inertia - node_.inertia)) ;

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
  
Fusion_ward1d::Fusion_ward1d(node_ward1d *node1_, node_ward1d *node2_) 
  : node1(node1_),  node2(node2_) 
  {
    double inertia_c1c2 = node1->x2plus + node2->x2plus - pow(node1->xplus + node2->xplus, 2) / (node1->weight + node2->weight);
    lambda = std::max(inertia_c1c2 - node1->inertia - node2->inertia, 0.0);
  }
;
