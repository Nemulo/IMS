#include <vector>
#include "cell.hpp"

#ifndef GRID
#define GRID

void create_grid(std::vector<std::vector<Cell>> &grid,double,double,double,double);

double distance_from_center(int,int,int,int);

void create_necrotic(std::vector<std::vector<Cell>> &grid,std::vector<std::vector<Cell>> &Newgrid);

void update_all_cells(std::vector<std::vector<Cell>> &grid,std::vector<std::vector<Cell>> &Newgrid);

void divide_proliferate(std::vector<std::vector<Cell>> &grid,std::vector<std::vector<Cell>> &NewGrid);

double calculate_average_tumor_radius(std::vector<std::vector<Cell>> &grid);

void create_initial_tumor(std::vector<std::vector<Cell>> &grid);

double calculate_average_necrotic_radius(std::vector<std::vector<Cell>> &grid);

std::vector<int> check_neighbours(std::vector<std::vector<Cell>> &grid,int,int,cell_type);

std::vector<int> find_closest_rim_cell(std::vector<std::vector<Cell>> &grid,int,int,cell_type);

void increment_time(std::vector<std::vector<Cell>> &grid, std::vector<std::vector<Cell>> &Newgrid);

bool print_info_tumor(std::vector<std::vector<Cell>> &grid,std::vector<std::vector<Cell>> &Newgrid,int);

#endif