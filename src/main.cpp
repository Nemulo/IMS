#include "cell.hpp"
#include <unistd.h>
#include <string>
#include <GL/glut.h>
#include <vector>
#include "grid.hpp"

#define MAX_GRID_SIZE 1000

std::vector<std::vector<Cell>> grid,NewGrid;



void print_help()
{
    std::string help = "Usage:\n\t./tumor -h\t-> show help\n./tumor -a <a> -b <b> -R <Rmax> -p <p0> [-g]\nRequired argumets:\n\t-a <a> -> float value (delimiter '.'): Base necrotic thickness, controlled by nutritional needs\n-b <b> float value (delimiter'.'): Base proliferative thickness\n-g -> Optional argument -> show automaton in gui\n";
    std::cout<<help<<std::endl;
}

// callback for GL
void display()
{
    glClear(GL_COLOR_BUFFER_BIT);
    glLoadIdentity();
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(0,1000,0,1000);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    glViewport(0, 0, MAX_GRID_SIZE, MAX_GRID_SIZE);


    for (const auto& row: grid)
    {
        for (const auto& cell: row)
        {
            switch (cell.ctype)
            {
                case cell_type::necrotic:
                    glColor3f(0.0,0.0,0.0);
                    break;

                case cell_type::proliferative:
                    glColor3f(1.0,0.0,0.0);
                    break;
                
                case cell_type::non_proliferative:
                    glColor3f(1.0,1.0,0.0);
                    break;

                case cell_type::healthy:
                    glColor3f(1.0,1.0,1.0);
                    break;
            }

            glBegin(GL_POINTS);
                glVertex2f(0.0f+cell.x_pos, 0.0f+cell.y_pos);
            glEnd();
        }
    }
    increment_time(grid,NewGrid);

    grid = NewGrid;

    glutSwapBuffers();
    glutPostRedisplay();
}

int main(int argc,char **argv)
{
    //-----------------------------------
    //              Argparse
    //-----------------------------------
    int opt;
    float a_factor,b_factor,probability, Rmax = 0;
    bool gui,a,b,p,r = false;
    bool err = false;
    while ((opt = getopt(argc,argv,"a:b:p:R:hg")) != -1)
    {
        switch (opt)
        {
            case 'a':
                a = true;
                if (optarg!=nullptr)
                {
                    std::string a_str = optarg;
                    try
                    {
                        a_factor = stof(a_str);
                    }
                    catch(const std::invalid_argument& ia)
                    {
                        std::cerr <<"Invalid argument 'a': "<< ia.what() << '\n';
                        err = true;
                    }
                    catch (const std::out_of_range& oor)
                    {
                        std::cerr <<"Out of range 'a':"<< oor.what()<<std::endl;
                        err = true;
                    }
                }
                break;
            case 'b':
                b = true;
                if (optarg!=nullptr)
                {
                    std::string b_str = optarg;
                    try
                    {
                        b_factor = stof(b_str);
                    }
                    catch(const std::invalid_argument& ia)
                    {
                        std::cerr <<"Invalid argument 'b': "<< ia.what() << '\n';
                        err = true;
                    }
                    catch (const std::out_of_range& oor)
                    {
                        std::cerr <<"Out of range 'b':"<< oor.what()<<std::endl;
                        err = true;
                    }
                }
                break;
            case 'R':
                r = true;
                if (optarg!=nullptr)
                {
                    std::string r_str = optarg;
                    try
                    {
                        Rmax = stof(r_str);
                        if (Rmax>50.0)
                        {
                            std::cerr<<"Rmax is limited to 50.0. Please rerun the program with correct value"<<std::endl;
                            err=true;
                        }
                    }
                    catch(const std::invalid_argument& ia)
                    {
                        std::cerr <<"Invalid argument 'R': "<< ia.what() << '\n';
                        err = true;
                    }
                    catch (const std::out_of_range& oor)
                    {
                        std::cerr <<"Out of range 'R':"<< oor.what()<<std::endl;
                        err = true;
                    }
                }
                break;
            case 'p':
                p = true;
                if (optarg!=nullptr)
                {
                    std::string p_str = optarg;
                    try
                    {
                        probability = stof(p_str);
                        if (probability>1.0 || probability<0)
                        {
                            std::cerr<<"Invalid probability arg, specify number 0-1"<<std::endl;
                            err=true;
                        }
                    }
                    catch(const std::invalid_argument& ia)
                    {
                        std::cerr <<"Invalid argument 'p': "<< ia.what() << '\n';
                        err = true;
                    }
                    catch (const std::out_of_range& oor)
                    {
                        std::cerr <<"Out of range 'p':"<< oor.what()<<std::endl;
                        err = true;
                    }
                    
                }
                break;
            case 'h':
                if (argc>2)
                {
                    std::cerr<<"Invalid argument combination"<<std::endl;
                    print_help();
                    return 1;
                }
                else
                {
                    print_help();
                    return 0;
                }
                break;
            case 'g':
                gui=true;
                break;
            default:
                print_help();
                return 1;
                break;
        }
    }

    if (err)
    {
        std::cerr<<"There was an error parsing arguments, the program will now exit"<<std::endl;
        return 1;
    }
    if (!a || !b || !p || !r)
    {
        std::cerr << "Missing required argument/s, please rerun the program with all required arguments"<<std::endl;
        print_help();
        return 1;
    }
    //----------------------------------------

    create_grid(grid, a_factor,b_factor,Rmax*10.0,probability);
    create_grid(NewGrid, a_factor,b_factor,Rmax*10.0,probability);
    create_initial_tumor(grid);
    NewGrid = grid;
    // visualizing:
    if (gui)
    {
        glutInit(&argc,argv);
        glutInitDisplayMode(GLUT_RGB);
        glutInitWindowPosition(0,0);
        glutInitWindowSize(1000,1000);

        glutCreateWindow("Brain tumor growth simulation");

        glutDisplayFunc(display);
        glClearColor(1.0,1.0,1.0,1.0);
        glutMainLoop();

    }
    if (!gui)
    {
        int day_count = 0;
        bool patient_alive = true;
        while (patient_alive)
        {
            patient_alive = print_info_tumor(grid,NewGrid,day_count);
            grid = NewGrid;
            day_count++;
        }
    }
    return 0;
}
