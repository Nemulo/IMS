#ifndef CELL_H
#define CELL_H

#include <iostream>
#include <stdlib.h>

enum struct cell_type
{
  healthy,
  non_proliferative,
  proliferative,
  necrotic
};


class Cell {
    private:
    /* data */
    public:

      int x_pos = 0;
      int y_pos = 0;
      double avgRt = 0;
      cell_type ctype;
      double a_arg = 0.0;
      double b_arg = 0.0;
      double Rmax = 0.0;
      double division_prob = 0.0;

      Cell(int,int, cell_type,double,double,double,double,double);

      cell_type ret_cell_type();
        
  };



#endif