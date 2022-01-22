#include "grid.hpp"
#include <math.h>
#include <iostream>
#include <time.h>
#include <unistd.h>
// initial radius of proliferative cells
#define INIT_RADIUS 30
#define CENTER 500

/*
    Function to search cell neighbours via Moore's neighbourhood
    
    params:
        grid -> 2D vector of cells to search in
        x,y -> coords of cell to be looked around
        ctype_searched -> type of cell that is current cell
    
    returns:
        vector of coords of the cell that is not the cell_type of ctype_searched
        or [-1,-1] if no such cell is found

*/
std::vector<int> check_neighbours(std::vector<std::vector<Cell>> &grid,int x,int y,cell_type ctype_searched)
{
    std::vector<int> searched_cell_coord;
    //resolve border cells
    if (x==0)
    {
        // left line
        if (y == 0)
        {
            //bottom left
            if (grid[x][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            searched_cell_coord.push_back(-1);
            searched_cell_coord.push_back(-1);
            return searched_cell_coord;
        }
        else if (y==999)
        {
            //top left
            if (grid[x][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            searched_cell_coord.push_back(-1);
            searched_cell_coord.push_back(-1);
            return searched_cell_coord;
        }
        else
        {
            //just left line
            if (grid[x][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            searched_cell_coord.push_back(-1);
            searched_cell_coord.push_back(-1);
            return searched_cell_coord;
        }
    }
    if (x==999)
    {
        //right line
        if (y == 0)
        {
            //bottom right
            if (grid[x-1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x-1][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            searched_cell_coord.push_back(-1);
            searched_cell_coord.push_back(-1);
            return searched_cell_coord;
        }
        else if (y==999)
        {
            // top right
            if (grid[x-1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x-1][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            searched_cell_coord.push_back(-1);
            searched_cell_coord.push_back(-1);
            return searched_cell_coord;
        }
        else
        {
            //just right line
            if (grid[x][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x-1][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x-1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x-1][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            searched_cell_coord.push_back(-1);
            searched_cell_coord.push_back(-1);
            return searched_cell_coord;
        }
    }
    if (y==0)
    {
        //just bottom line
            if (grid[x-1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x-1][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y+1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            searched_cell_coord.push_back(-1);
            searched_cell_coord.push_back(-1);
            return searched_cell_coord;
    }
    if (y==999)
    {
        //just top line
            if (grid[x-1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x-1][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y-1].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            if (grid[x+1][y].ctype != ctype_searched)
            {
                searched_cell_coord.push_back(x);
                searched_cell_coord.push_back(y);
                return searched_cell_coord;
            }
            searched_cell_coord.push_back(-1);
            searched_cell_coord.push_back(-1);
            return searched_cell_coord;
    }
    // non-border cells
    if (grid[x-1][y].ctype != ctype_searched)
    {
        searched_cell_coord.push_back(x);
        searched_cell_coord.push_back(y);
        return searched_cell_coord;
    }
    if (grid[x-1][y+1].ctype != ctype_searched)
    {
        searched_cell_coord.push_back(x);
        searched_cell_coord.push_back(y);
        return searched_cell_coord;
    }
    if (grid[x][y+1].ctype != ctype_searched)
    {
        searched_cell_coord.push_back(x);
        searched_cell_coord.push_back(y);
        return searched_cell_coord;
    }
    if (grid[x+1][y+1].ctype != ctype_searched)
    {
        searched_cell_coord.push_back(x);
        searched_cell_coord.push_back(y);
        return searched_cell_coord;
    }
    if (grid[x+1][y].ctype != ctype_searched)
    {
        searched_cell_coord.push_back(x);
        searched_cell_coord.push_back(y);
        return searched_cell_coord;
    }
    if (grid[x+1][y-1].ctype != ctype_searched)
    {
        searched_cell_coord.push_back(x);
        searched_cell_coord.push_back(y);
        return searched_cell_coord;
    }
    if (grid[x][y-1].ctype != ctype_searched)
    {
        searched_cell_coord.push_back(x);
        searched_cell_coord.push_back(y);
        return searched_cell_coord;
    }
    if (grid[x-1][y-1].ctype != ctype_searched)
    {
        searched_cell_coord.push_back(x);
        searched_cell_coord.push_back(y);
        return searched_cell_coord;
    }
    searched_cell_coord.push_back(-1);
    searched_cell_coord.push_back(-1);
    return searched_cell_coord;
}


/*
    Function to create initial tumor on grid of size INIT_RADIUS

    params:
        grid -> grid that will have the initial tumor inserted

*/
void create_initial_tumor(std::vector<std::vector<Cell>> &grid)
{
    for (int i = 0 ; i<1000 ; i++)
    {
        for (int j = 0 ; j<1000 ; j++)
        {
            if (distance_from_center(i,j,CENTER,CENTER) <= INIT_RADIUS)
            {
                grid[i][j].ctype = cell_type::proliferative;
                grid[i][j].avgRt = INIT_RADIUS;
            }
            else
            {
                grid[i][j].ctype = cell_type::healthy;
                grid[i][j].avgRt = INIT_RADIUS;
            }
        }
    }
}

/*
    Function to create clean grid of cells

    params:
        grid -> cells will be inserted into this 2D vector
        a,b,rmax,p0 -> every cell will hold these constants from program params

*/
void create_grid(std::vector<std::vector<Cell>> &grid,double a,double b,double rmax,double p0)
{
    srand(time(NULL));
    // for every row, for every col in row, create new cell with x of iteration,y of iteration, use distance from center
    //  if distance from center <50, set its status as proliferative, otherwise healthy
    for (int i = 0 ; i<1000 ; i++)
    {
        std::vector<Cell> cell_row {};
        for (int j = 0 ; j<1000 ; j++)
        {
            cell_row.push_back(Cell(i,j,cell_type::healthy,0.0,a,b,rmax,p0));
        }
        grid.push_back(cell_row);
    }
}


/*
    Function to calculate distance of current cell from center cell (or any cell)
    Uses pythagorean theorem.

    Params:
        x_curr,y_curr -> coords of current cell
        x_cent,y_cent -> corrds of cell to calculate distance from

    Returns:
        distance of current cell to center cell

*/
double distance_from_center(int x_curr,int y_curr, int x_cent, int y_cent)
{
    //pythagorean theorem
    double a = x_curr - x_cent;
    double b = y_curr - y_cent;

    a = pow(a,2.0);
    b = pow(b,2.0);

    return sqrt(a+b);
}


/* 
    Function To search all non-proliferate cells, apply rules on them and change them to necrotic if needed
    Function also copies all existing necrotic cells to new grid

    Params:
        grid -> current grid of cells
        Newgrid -> grid of cells in next time iteration

*/
void create_necrotic(std::vector<std::vector<Cell>> &grid,std::vector<std::vector<Cell>> &Newgrid)
{
    //save non-prolif cells
    std::vector<Cell>non_prolif_cells;
    for (int i = 0 ; i < 1000 ; i++)
    {
        for (int j = 0 ; j < 1000 ; j ++)
        {
            if (grid[i][j].ctype == cell_type::non_proliferative)
            {
                //save cell
                non_prolif_cells.push_back(grid[i][j]);
            }
        }
    }

    //process non-prolif cells
    for (auto cell : non_prolif_cells)
    {
        //find closest proliferate rim cell
        std::vector<int> closest_prolif_rim_cell = find_closest_rim_cell(grid,cell.x_pos,cell.y_pos,cell_type::proliferative);
        //calculate max range to stay non prolif
        double delta_n = (cell.a_arg*pow(cell.avgRt,2.0/3.0));
        if (distance_from_center(cell.x_pos,cell.y_pos,closest_prolif_rim_cell[0],closest_prolif_rim_cell[1])>delta_n)
        {
            Newgrid[cell.x_pos][cell.y_pos] = cell;
            Newgrid[cell.x_pos][cell.y_pos].ctype = cell_type::necrotic;
        }
        else
        {
            Newgrid[cell.x_pos][cell.y_pos] = cell;
        }
    }

    //save all necrotic
    std::vector<Cell>necro;
    for (int i = 0 ; i < 1000 ; i++)
    {
        for (int j = 0 ; j < 1000 ; j ++)
        {
            if (grid[i][j].ctype == cell_type::necrotic)
            {
                //save cell
                necro.push_back(grid[i][j]);
            }
        }
    }
    for (auto cell : necro)
    {
        Newgrid[cell.x_pos][cell.y_pos] = cell;
    }
}

/*
    Function to calculate average tumor radius, for more info see documentation.pdf

    Params:
        grid -> current grid of cells 

    Returns:
        average tumor radius in cells

*/
double calculate_average_tumor_radius(std::vector<std::vector<Cell>> &grid)
{
    std::vector<double> prolif_cells_radius;
    for (int i=0 ; i<1000 ; i++)
    {
        for (int j=0 ; j<1000 ; j++)
        {
            if (grid[i][j].ctype == cell_type::proliferative)
            {
                std::vector<int> coords = check_neighbours(grid,i,j,cell_type::proliferative);
                if (coords[0] == -1)
                {
                    continue;
                }
                else
                {
                    double dist = distance_from_center(coords[0],coords[1],CENTER,CENTER);
                    prolif_cells_radius.push_back(dist);
                }
            }
        }
    }
    double Rt = 0.0;
    for (auto value : prolif_cells_radius)
    {
        Rt = Rt + value;
    }
    if (prolif_cells_radius.size()!=0)
    {   
        Rt = Rt/prolif_cells_radius.size();
        return Rt;
    }
    return 0.0;
}

/*
    Same function as calculate_average_tumor radius, only this one focuses on necrotic core
    See function above for more details

*/
double calculate_average_necrotic_radius(std::vector<std::vector<Cell>> &grid)
{
    std::vector<double> necro_cells_radius;
    for (int i=0 ; i<1000 ; i++)
    {
        for (int j=0 ; j<1000 ; j++)
        {
            if (grid[i][j].ctype == cell_type::necrotic)
            {
                std::vector<int> coords = check_neighbours(grid,i,j,cell_type::necrotic);
                if (coords[0] == -1)
                {
                    continue;
                }
                else
                {
                    double dist = distance_from_center(coords[0],coords[1],CENTER,CENTER);
                    necro_cells_radius.push_back(dist);
                }
            }
        }
    }
    double Rt = 0.0;
    for (auto value : necro_cells_radius)
    {
        Rt = Rt + value;
    }
    if (necro_cells_radius.size()!= 0 )
    {
        Rt = Rt/necro_cells_radius.size();
        return Rt; 
    }
    return 0.0;
}

/*
    Function for applying all rules on all cells in current grid
    
    params:
        grid -> current grid of cells
        Newgrid -> grid of cells in next iteration
*/
void update_all_cells(std::vector<std::vector<Cell>> &grid,std::vector<std::vector<Cell>> &Newgrid)
{
    
    divide_proliferate(grid,Newgrid);
    
    create_necrotic(grid,Newgrid);
}

/*
    Function to create new proliferative cells in grid/change proliferative to nonproliferative cells according to existing rules(see documentation.pdf)

    params:
        grid -> grid of current cells
        Newgrid -> grid of cells in next iteration

*/
void divide_proliferate(std::vector<std::vector<Cell>> &grid,std::vector<std::vector<Cell>> &NewGrid)
{
    //save prolif cells
    std::vector<Cell>prolif_cells;
    for (int i = 0 ; i < 1000 ; i++)
    {
        for (int j = 0 ; j < 1000 ; j ++)
        {
            if (grid[i][j].ctype == cell_type::proliferative)
            {
                //save cell
                prolif_cells.push_back(grid[i][j]);
            }
        }
    }
    
    //process prolif cells
    for (auto cell:prolif_cells)
    {
        double Pd = 0.0;
        Pd = cell.division_prob*(1-(distance_from_center(cell.x_pos,cell.y_pos,CENTER,CENTER)/cell.Rmax));
        Pd = round(Pd * 1000.0);
        int div = rand() % 1000 + 1;
        if (div<Pd)
        {
            //save actuall cell to new grid
            NewGrid[cell.x_pos][cell.y_pos] = cell;
            
            
            // divide cell
            //find division radius max
            double delta_p = cell.b_arg*pow(cell.avgRt,2.0/3.0);
            // search always closest radius
            int radius = 0;
            std::vector<Cell> healthy_cells;
            while(radius<=delta_p+1.0)
            {
                
                for (int i = cell.x_pos-radius; i <= cell.x_pos+radius ; i++)
                {
                    for (int j = cell.y_pos-radius; j<= cell.y_pos+radius ; j++)
                    {
                        if(distance_from_center(grid[i][j].x_pos,grid[i][j].y_pos,cell.x_pos,cell.y_pos)<delta_p)
                        {
                            
                            if (grid[i][j].ctype == cell_type::healthy)
                            {
                                healthy_cells.push_back(grid[i][j]);
                            }
                        }
                    }
                }
                
                if (healthy_cells.size()!=0)
                {
                    int div_cell = rand() % healthy_cells.size();
                    NewGrid[healthy_cells[div_cell].x_pos][healthy_cells[div_cell].y_pos] = grid[healthy_cells[div_cell].x_pos][healthy_cells[div_cell].y_pos];
                    NewGrid[healthy_cells[div_cell].x_pos][healthy_cells[div_cell].y_pos].ctype = cell_type::proliferative;
                    break;
                }
                radius++;
            }
            
            if (healthy_cells.size() == 0)
            {
                NewGrid[cell.x_pos][cell.y_pos].ctype = cell_type::non_proliferative;
            }
            
        }
    }

}

/*
    Function to find closest cell of type ctype_searched (used for applying some rule(s), see documentation.pdf)

    params:
        grid -> grid of current cells
        x,y -> coords of current cell
        ctype_searched -> cell_type of cell searched

    returns:
        coords of closest cell of ctype_searched

*/
std::vector<int> find_closest_rim_cell(std::vector<std::vector<Cell>> &grid,int x,int y,cell_type ctype_searched)
{

    for (int i = 0 ; ; i++)
    {
       for (int cell_rows = x-i ; cell_rows <= x+i ; cell_rows++)
       {
           for (int cell_cols = y-i ; cell_cols <= y+i ; cell_cols++)
           {
               int dist = distance_from_center(cell_rows,cell_cols,x,y);
               if (dist == i)
               {
                   if ((cell_rows < 0 || cell_rows >= 1000) || (cell_cols<0 || cell_cols >= 1000))
                   {
                       continue;
                   }
                   else
                   {
                       if (grid[cell_rows][cell_cols].ctype == ctype_searched)
                       {
                           std::vector<int> rim_cell;
                           rim_cell.push_back(cell_rows);
                           rim_cell.push_back(cell_cols);
                           return rim_cell;
                       }
                   }
               }
           }
       }

    }
}

/*
    Function to increment discrete type (proceed to next iteration of evolution of the tumor)
    with every iteration new average tumor size is inserted into every cell

    params:
        grid -> grid of current cells
        Newgrid -> grid of cells in next iteration


*/
void increment_time(std::vector<std::vector<Cell>> &grid, std::vector<std::vector<Cell>> &Newgrid)
{
    double Rt = calculate_average_tumor_radius(grid);
    for (int i = 0 ; i < 1000 ; i++)
    {
        for (int j = 0 ; j < 1000 ; j++)
        {
            grid[i][j].avgRt = Rt;
        }
    }
    update_all_cells(grid,Newgrid);
}

/*
    Function to print text info about tumor if GUI option was not selected
    uses increment time to run simulation similiarly to gui option

    params:
        grid -> current grid of cells
        NewGrid -> grid of cells in next iteration
        day_cnt -> count of days

    Returns:
        true if patient is still alive, false otherwise


*/
bool print_info_tumor(std::vector<std::vector<Cell>> &grid,std::vector<std::vector<Cell>> &NewGrid, int day_cnt)
{
    increment_time(grid,NewGrid);
    double volume = (4.0/3.0)*3.14*pow(grid[0][0].avgRt/10.0,3.0);
    if ((day_cnt % 25) == 0)
    {
        std::cout<<"\t#\t#\t# Day "<<day_cnt<<std::endl;
        double necro_rad = calculate_average_necrotic_radius(grid);
        std::cout<<"Average Tumor Size: "<<grid[0][0].avgRt/10.0<<" mm\tAverge necrotic core radius :"<<necro_rad/10.0<<" mm"<<std::endl;
        int tumor_cells = 0;
        int necro_cells = 0;
        for (int i = 0; i<1000 ; i++)
        {
            for (int j = 0 ; j < 1000 ; j++)
            {
                if (grid[i][j].ctype == cell_type::non_proliferative || grid[i][j].ctype == cell_type::proliferative)
                {
                    tumor_cells++;
                }
                if(grid[i][j].ctype == cell_type::necrotic)
                {
                    necro_cells++;
                }
            }
        }
        std::cout<<"Cell number: \n\tnecrotic: "<<necro_cells<<"\ttumor: "<<tumor_cells<<std::endl;
        std::cout<<"Tumor volume: "<<volume/1000.0<<" cubic cm"<<std::endl;
    }
    if (volume>65000 && day_cnt>25)
    {
        std::cout<<"\t#\t#\t# Day "<<day_cnt<<std::endl;
        double necro_rad = calculate_average_necrotic_radius(grid);
        std::cout<<"Average Tumor Size: "<<grid[0][0].avgRt/10.0<<" mm\tAverge necrotic core radius :"<<necro_rad/10.0<<std::endl;
        int tumor_cells = 0;
        int necro_cells = 0;
        for (int i = 0; i<1000 ; i++)
        {
            for (int j = 0 ; j < 1000 ; j++)
            {
                if (grid[i][j].ctype == cell_type::non_proliferative || grid[i][j].ctype == cell_type::proliferative)
                {
                    tumor_cells++;
                }
                if(grid[i][j].ctype == cell_type::necrotic)
                {
                    necro_cells++;
                }
            }
        }
        std::cout<<"Cell number: \n\tnecrotic: "<<necro_cells<<"\ttumor: "<<tumor_cells<<std::endl;
        std::cout<<"Tumor volume: "<<volume/1000.0<<" cubic cm"<<std::endl;
        std::cout<<"Tumor has reached volume of 65 cubic cm, patient is dead\n";
        return false;
    }
    return true;
}