#include "cell.hpp"

/*
    Object constructor
        Creates new cell and assigns all neccesary parameters as object variables

    params:
        x,y -> position on grid
        celltype -> cell_type of current cell
        avgSize -> Rt of tumor
        a,b,max,p0 -> constants from program params input
*/
Cell::Cell(int x, int y, cell_type celltype,double avgSize,double a,double b,double max,double p0)
{
    this->x_pos = x;
    this->y_pos = y;
    this->ctype = celltype;
    this->avgRt = avgSize;
    this->a_arg = a;
    this->b_arg = b;
    this->Rmax = max;
    this->division_prob = p0;
}

/*
    Unused function
        serves for debbugging and checking cell types

*/
cell_type Cell::ret_cell_type()
{
    return this->ctype;
}