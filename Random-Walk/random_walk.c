/* By Aleksandar Mitic */
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <time.h>
#include <math.h>

#define NPROC0 4
#define NPROC1 4
#define min(X, Y) (((X) <= (Y)) ? (X) : (Y))

#define L0 40
#define L1 40
#define L2 40
#define VOLUME L0 *L1 *L2
#define iL0 L0 / NPROC0
#define iL1 L1 / NPROC1

#define UP 0
#define DOWN 1
#define LEFT 2
#define RIGHT 3

typedef struct GRID_INFO_T
{
    int p;
    MPI_Comm comm;
    MPI_Comm my_row;
    MPI_Comm my_col;
    int q;
    int boundary;
    int grid_rank;
    int work_rank;
    int cpr[2];
    int npr[4];
    int face0;
    int face1;
    int gl_coord[3];
    int starting_global_coordinates[3];
    int ipt[iL0 * iL1 * L2];
    int iup[iL0 * iL1 * L2][2];
    int idn[iL0 * iL1 * L2][2];
    int *map;
    int cur_rank;

} GRID_INFO_T;

void Setup_grid(GRID_INFO_T *grid)
{
    grid->face0 = (L2 * iL0);
    grid->face1 = (L2 * iL1);

    grid->boundary = 2 * (grid->face0 + grid->face1);
    grid->map = (int *)malloc(grid->boundary * sizeof(int));

    int periods[2] = {1, 1}; /* Set the lattice to be periodic */
    int free_coords[2];

    MPI_Comm_size(MPI_COMM_WORLD, &(grid->p));
    grid->q = sqrt((NPROC0 * NPROC1));

    if (grid->p < 2)
    {
        printf("Number of %d processors not enough\n", grid->p);
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    /* Create a 2D grid structure */
    int dims[2] = {grid->q, grid->q};

    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &(grid->comm));
    /* MPI_Comm_rank(grid->comm, &(grid->my_rank)); */

    /* Get the rank and coordinates of the process in the lattice */
    MPI_Comm_rank(grid->comm, &(grid->grid_rank));
    MPI_Cart_coords(grid->comm, grid->grid_rank, 2, grid->cpr);

    /* Get the ranks of the neighboring processes */
    MPI_Cart_shift(grid->comm, 0, 1, &(grid->npr[DOWN]), &(grid->npr[UP]));
    MPI_Cart_shift(grid->comm, 1, 1, &(grid->npr[LEFT]), &(grid->npr[RIGHT]));
    /* MPI_Comm_free(&grid->comm); */
}
void define_ipt(GRID_INFO_T *grid) /* Function to calculate ipt */
{
    int vol;
    vol = iL0 * iL1 * L2;
    int beg = 0;
    int mid = (iL0 * iL1 * L2) / 2;

    for (int i = 0; i < L2; i++)
    {
        for (int j = 0; j < iL1; j++)
        {
            for (int k = 0; k < iL0; k++)
            {
                if ((i + j + k) % 2 == 0)
                {
                    grid->ipt[k + j * iL0 + i * iL0 * iL1] = beg;
                    beg++;
                }
                else
                {
                    grid->ipt[k + j * iL0 + i * iL0 * iL1] = mid;
                    mid++;
                }
            }
        }
    }
}

void define_iup(GRID_INFO_T *grid)
{
    int vol;
    vol = iL0 * iL1 * L2;
    int bndry = 2 * (grid->face0 + grid->face1);
    int b = 0;
    int count[4] = {0, 0, 0, 0}; /* for updating every element in iup for all 4 bounday conditions */

    /* Setting up the index values for positive 0 and 1 direction for even and odd points */
    int *pos0_o, *pos0_e;
    int *pos1_o, *pos1_e;

    pos0_o = (int *)malloc(grid->face0 / 2 * sizeof(int));
    pos0_e = (int *)malloc(grid->face0 / 2 * sizeof(int));
    pos1_o = (int *)malloc(grid->face1 / 2 * sizeof(int));
    pos1_e = (int *)malloc(grid->face1 / 2 * sizeof(int));

    int cnt = 0;
    for (int i = vol + grid->face0 / 2; i < vol + grid->face0; i++)
    {
        pos0_e[cnt] = i;
        cnt++;
    }
    cnt = 0;
    for (int i = vol + grid->face0 + grid->face1 / 2; i < vol + bndry / 2; i++)
    {
        pos1_e[cnt] = i;
        cnt++;
    }
    cnt = 0;
    for (int i = vol + bndry / 2 + grid->face0 / 2; i < vol + bndry / 2 + grid->face0; i++)
    {
        pos0_o[cnt] = i;
        cnt++;
    }
    cnt = 0;
    for (int i = vol + bndry / 2 + grid->face0 + grid->face1 / 2; i < vol + bndry; i++)
    {
        pos1_o[cnt] = i;
        cnt++;
    }

    int sizet = sizeof(pos0_e) / sizeof(int);
    /* printf("size of poso_e is %d\n",sizet); */

    int ix;
    for (int k = 0; k < L2; k++)
    {
        for (int i = 0; i < iL0; i++)
        {
            for (int j = 0; j < iL1; j++)
            {
                ix = i + j * iL0 + k * iL0 * iL1; /* Current lexicographic index*/
                if (j == iL1 - 1 && i < iL0 - 1)  /* Last column*/
                {
                    if ((i + j + k) % 2 == 0)
                    {
                        grid->iup[ix][1] = pos1_e[count[0]];
                        grid->iup[ix][0] = ix + 1;
                        count[0]++;
                    }
                    else
                    {
                        grid->iup[ix][1] = pos1_o[count[1]];
                        grid->iup[ix][0] = ix + 1;
                        count[1]++;
                    }
                }
                else if (i == iL0 - 1 && j < iL1 - 1) /* Last row */
                {
                    if ((i + j + k) % 2 == 0)
                    {
                        grid->iup[ix][0] = pos0_e[count[2]];
                        grid->iup[ix][1] = ix + iL0;
                        count[2]++;
                    }
                    else
                    {
                        grid->iup[ix][0] = pos0_o[count[3]];
                        grid->iup[ix][1] = ix + iL0;
                        count[3]++;
                    }
                }
                else if (i == iL0 - 1 && j == iL0 - 1) /* top right corner element in the local lattice which has both iup values outside the local lattice */
                {
                    if ((i + j + k) % 2 == 0)
                    {
                        grid->iup[ix][0] = pos0_e[count[2]];
                        grid->iup[ix][1] = pos1_e[count[0]];
                        count[2]++;
                        count[0]++;
                    }
                    else
                    {
                        grid->iup[ix][0] = pos0_o[count[3]];
                        grid->iup[ix][1] = pos1_o[count[1]];
                        count[3]++;
                        count[1]++;
                    }
                }
                else /* values inside the lattice which have iup within the local lattice */
                {
                    grid->iup[ix][0] = ix + 1;
                    grid->iup[ix][1] = ix + iL0;
                }
            }
        }
    }
    free(pos0_o);
    free(pos0_e);
    free(pos1_e);
    free(pos1_o);
}

void define_idn(GRID_INFO_T *grid)
{
    int vol;
    vol = iL0 * iL1 * L2;
    int bndry = 2 * (grid->face0 + grid->face1);
    int b = 0;
    int count[4] = {0, 0, 0, 0}; /* for updating every element in iup for all 4 bounday conditions*/

    int *neg0_o, *neg0_e;
    int *neg1_o, *neg1_e;

    neg0_o = (int *)malloc(grid->face0 / 2 * sizeof(int));
    neg0_e = (int *)malloc(grid->face0 / 2 * sizeof(int));
    neg1_o = (int *)malloc(grid->face1 / 2 * sizeof(int));
    neg1_e = (int *)malloc(grid->face1 / 2 * sizeof(int));

    int cnt = 0;
    for (int i = vol; i < vol + grid->face0 / 2; i++)
    {
        neg0_e[cnt] = i;
        cnt++;
    }
    cnt = 0;
    for (int i = vol + grid->face0; i < vol + grid->face0 + grid->face1 / 2; i++)
    {
        neg1_e[cnt] = i;
        cnt++;
    }
    cnt = 0;
    for (int i = vol + bndry / 2; i < vol + bndry / 2 + grid->face0 / 2; i++)
    {
        neg0_o[cnt] = i;
        cnt++;
    }
    cnt = 0;
    for (int i = vol + bndry / 2 + grid->face0; i < vol + bndry / 2 + grid->face0 + grid->face1 / 2; i++)
    {
        neg1_o[cnt] = i;
        cnt++;
    }

    int ix;
    for (int k = 0; k < L2; k++)
    {
        for (int i = 0; i < iL0; i++)
        {
            for (int j = 0; j < iL1; j++)
            {
                ix = i + j * iL0 + k * iL0 * iL1; /* Current lexicographic index*/
                if (j == 0 && i > 0)              /* first column*/
                {
                    if ((i + j + k) % 2 == 0)
                    {
                        grid->idn[ix][1] = neg1_e[count[0]];
                        grid->idn[ix][0] = ix - 1;
                        count[0]++;
                    }
                    else
                    {
                        grid->idn[ix][1] = neg1_o[count[1]];
                        grid->idn[ix][0] = ix - 1;
                        count[1]++;
                    }
                }
                else if (i == 0 && j > 0) /* first row */
                {
                    if ((i + j + k) % 2 == 0)
                    {
                        grid->idn[ix][0] = neg0_e[count[2]];
                        grid->idn[ix][1] = ix - iL0;
                        count[2]++;
                    }
                    else
                    {
                        grid->idn[ix][0] = neg0_o[count[3]];
                        grid->idn[ix][1] = ix - iL0;
                        count[3]++;
                    }
                }
                else if (i == 0 && j == 0) /* bottom-left element in the local lattice which has both idn values outside the local lattice */
                {
                    if ((i + j + k) % 2 == 0)
                    {
                        grid->idn[ix][0] = neg0_e[count[2]];
                        grid->idn[ix][1] = neg1_e[count[0]];
                        count[2]++;
                        count[0]++;
                    }
                    else
                    {
                        grid->idn[ix][0] = neg0_o[count[3]];
                        grid->idn[ix][1] = neg1_o[count[1]];
                        count[3]++;
                        count[1]++;
                    }
                }
                else /* values inside the lattice which have iup within the local lattice */
                {
                    grid->idn[ix][0] = ix - 1;
                    grid->idn[ix][1] = ix - iL0;
                }
            }
        }
    }
    free(neg0_o);
    free(neg0_e);
    free(neg1_e);
    free(neg1_o);
}

void define_map(GRID_INFO_T *grid)
{
    int loc_vol = iL0 * iL1 * L2;
    int bndry = 2 * (grid->face0 + grid->face1); // 400 points
    int ix;
    for (int k = 0; k < L2; k++)
    {
        for (int i = 0; i < iL0; i++)
        {
            for (int j = 0; j < iL1; j++)
            {
                ix = i + j * iL0 + k * iL0 * iL1; /* Current lexicographic index*/

                /* Verifying iup and idn to map local indices in neighbouring process*/

                /* Calculates the local indices in neighbouring process which will be a copy of subsequent element(first element) in the same row*/
                if (grid->iup[ix][0] >= loc_vol)
                {
                    // grid->map[grid->iup[ix][0]-loc_vol] = ix - (iL0-1);
                    grid->map[grid->iup[ix][0] - loc_vol] = ix - (iL0 - 1);
                }

                /* Calculates the local indices in neighbouring process which will be a copy of subsequent element(first element) in the same column*/
                if (grid->iup[ix][1] >= loc_vol)
                {
                    grid->map[grid->iup[ix][1] - loc_vol] = ix - ((iL0 - 1) * iL1);
                }
                /* Similar calculation for idn*/
                if (grid->idn[ix][0] >= loc_vol)
                {
                    grid->map[grid->idn[ix][0] - loc_vol] = ix + (iL0 - 1);
                }

                if (grid->idn[ix][1] >= loc_vol)
                {
                    grid->map[grid->idn[ix][1] - loc_vol] = ix + ((iL0 - 1) * iL1);
                }
            }
        }
    }
}

void global_coordinates_cal(GRID_INFO_T *grid, int *position)
{
    int NPROC_mu[3] = {grid->q, grid->q, 1};
    int L_mu[3] = {L0 / NPROC0, L1 / NPROC1, L2};
    grid->gl_coord[0] = L_mu[0] * grid->cpr[0] + position[0];
    grid->gl_coord[1] = L_mu[1] * grid->cpr[1] + position[1];
    grid->gl_coord[2] = position[2];

    grid->starting_global_coordinates[0] = L_mu[0] * grid->cpr[0] + position[0];
    grid->starting_global_coordinates[1] = L_mu[1] * grid->cpr[1] + position[1];
    grid->starting_global_coordinates[2] = position[2];

    grid->work_rank = grid->cpr[1] * NPROC0 + grid->cpr[0];
}

void check_coordinates(GRID_INFO_T *grid)
{
    if (grid->gl_coord[0] == L0)
    {
        grid->gl_coord[0] = 0;
    }
    if (grid->gl_coord[0] < 0)
    {
        grid->gl_coord[0] = L0 - 1;
    }
    if (grid->gl_coord[1] == L1)
    {
        grid->gl_coord[1] = 0;
    }
    if (grid->gl_coord[1] < 0)
    {
        grid->gl_coord[1] = L1 - 1;
    }
    if (grid->gl_coord[2] == L2)
    {
        grid->gl_coord[2] = 0;
    }
    if (grid->gl_coord[2] < 0)
    {
        grid->gl_coord[2] = L2 - 1;
    }
}

void distance_calculation(int *starting_position, int *new_position, int *NPROC_mu, int *L_mu, int mu, int i, int j)
{
    int N_meas = pow(10, 4);
    int T = 100;
    double dist1, dist2;
    double d_it[N_meas][T];
    if (mu == 0)
    {
        L_mu[0] = L0;
        L_mu[1] = 0;
        L_mu[2] = 0;
        NPROC_mu[0] = 1;
        NPROC_mu[1] = 0;
        NPROC_mu[2] = 0;
    }
    else if (mu == 1)
    {
        L_mu[0] = 0;
        L_mu[1] = L1;
        L_mu[2] = 0;
        NPROC_mu[0] = 0;
        NPROC_mu[1] = 1;
        NPROC_mu[2] = 0;
    }
    else
    {
        L_mu[0] = 0;
        L_mu[1] = 0;
        L_mu[2] = L2;
        NPROC_mu[0] = 0;
        NPROC_mu[1] = 0;
        NPROC_mu[2] = 1;
    }
    dist1 = pow(starting_position[0] - new_position[0], 2) + pow(starting_position[1] - new_position[1], 2) + pow(starting_position[2] - new_position[2], 2);
    dist2 = pow(NPROC_mu[0] * L_mu[0] - new_position[0] + starting_position[0], 2) + pow(NPROC_mu[1] * L_mu[1] - new_position[1] + starting_position[1], 2) + pow(NPROC_mu[2] * L_mu[2] - new_position[2] + starting_position[2], 2);
    d_it[j][i] = sqrt(min(dist1, dist2));
    printf("%f\n",d_it[j][i]);
    // //open a file for writing in append mode
    //FILE *fp = fopen("distance.txt", "a");
    // fprintf(fp, "%f\n", d_it[j][i]);
    // // close the file
    // fclose(fp);
}

void Random_Walk(GRID_INFO_T *grid)
{
    int direction;
    int mu;          /* new */
    int up_down_dir; /* new */
    int nl;          /* new */
    int shift_dir = 0;
    int i, t;
    int positon_start[3] = {0, 0, 0};
    int N_meas = pow(10, 4);
    int T = 100;
    int NPROC_mu[3] = {grid->q, grid->q, 1};
    int L_mu[3] = {L0 / NPROC0, L1 / NPROC1, L2};
    int loc_vol = iL0 * iL1 * L2;
    int rand_rank;
    MPI_Status status;

    /* Perform a random walk */
    for (i = 0; i < N_meas; i++)
    {
        if (grid->grid_rank==0)
        {
            rand_rank= (rand() % (grid->p));
        }
        MPI_Bcast(&rand_rank, 1, MPI_INT, 0, MPI_COMM_WORLD);

        if (grid->grid_rank == rand_rank)
        {
            /* Make initial guess for the starting position */
            positon_start[0] = (rand() % iL0);
            positon_start[1] = (rand() % iL1);
            positon_start[2] = (rand() % L2);
            nl = positon_start[0] + positon_start[1] * iL0 + positon_start[2] * iL0 * iL1;
            global_coordinates_cal(grid, positon_start);
        }
        MPI_Bcast(&(grid->work_rank), 1, MPI_INT, rand_rank, MPI_COMM_WORLD);
        MPI_Bcast(&(grid->starting_global_coordinates), 3, MPI_INT, rand_rank, MPI_COMM_WORLD);
        MPI_Bcast(&nl, 1, MPI_INT, rand_rank, MPI_COMM_WORLD);

        grid->gl_coord[0]= grid->starting_global_coordinates[0];
        grid->gl_coord[1]= grid->starting_global_coordinates[1];
        grid->gl_coord[2]= grid->starting_global_coordinates[2];

        grid->cur_rank = grid->work_rank;
        for (t = 0; t < T; t++)
        {
            if (grid->grid_rank == grid->work_rank)
            {
                /* Choose a random direction for the walk */
                direction = (rand() % 6);
                switch (direction) /* Moving +1 or -1 in one of the 6 directions */
                {
                case 0:
                    grid->gl_coord[0]++;
                    mu = 0;
                    shift_dir = 3; /* MPI_Cart_Shift */
                    break;
                case 1:
                    grid->gl_coord[0]--;
                    mu = 0;
                    shift_dir = 2; /* MPI_Cart_Shift */
                    break;
                case 2:
                    grid->gl_coord[1]++;
                    mu = 1;
                    shift_dir = 0; /* MPI_Cart_Shift */
                    break;
                case 3:
                    grid->gl_coord[1]--;
                    mu = 1;
                    shift_dir = 1; /* MPI_Cart_Shift */
                    break;
                case 4:
                    grid->gl_coord[2]++;
                    mu = 2; /* 3rd dimension */
                    up_down_dir = 1;
                    break;
                case 5:
                    grid->gl_coord[2]--;
                    mu = 2; /* 3rd dimension */
                    up_down_dir = 0;
                    break;
                default:
                    break;
                }
                check_coordinates(grid);
                if (mu != 2)
                {
                    if (shift_dir == 3) /* right */
                    {
                        if (grid->iup[nl][mu] < loc_vol)
                        {
                            nl = nl++;
                        }
                        else
                        {
                            nl = grid->map[grid->iup[nl][mu] - loc_vol];
                            grid->work_rank = grid->npr[3];
                        }
                    }
                    else if (shift_dir == 0) /* upward */
                    {
                        if (grid->iup[nl][mu] < loc_vol)
                        {
                            nl = nl + iL1;
                        }
                        else
                        {
                            nl = grid->map[grid->iup[nl][mu] - loc_vol];
                            grid->work_rank = grid->npr[0];
                        }
                    }
                    else if (shift_dir == 2) /* leftward */
                    {
                        if (grid->idn[nl][mu] < loc_vol)
                        {
                            nl = nl--;
                        }
                        else
                        {
                            nl = grid->map[grid->idn[nl][mu] - loc_vol];
                            grid->work_rank = grid->npr[2];
                        }
                    }
                    else if (shift_dir == 1) /* downward */
                    {
                        if (grid->idn[nl][mu] < loc_vol)
                        {
                            nl = nl - iL1;
                        }
                        else
                        {
                            nl = grid->map[grid->idn[nl][mu] - loc_vol];
                            grid->work_rank = grid->npr[1];
                        }
                    }
                }
                else
                {
                    if (up_down_dir == 1) /* top */
                    {
                        if (nl + iL0 * iL1 >= loc_vol)
                        {
                            nl = nl - (loc_vol - iL0 * iL1);
                        }
                        else
                        {
                            nl = nl + iL0 * iL1;
                        }
                    }
                    else /* bottom */
                    {
                        if (nl >= iL0 * iL1)
                        {
                            nl = nl - iL0 * iL1;
                        }
                        else
                        {
                            nl = nl + (loc_vol - iL1 * iL0);
                        }
                    }
                }
                distance_calculation(grid->starting_global_coordinates, grid->gl_coord, NPROC_mu, L_mu, mu, i, t);
                MPI_Barrier(MPI_COMM_WORLD);
            }
            else
            {
                MPI_Barrier(MPI_COMM_WORLD);
            }
            MPI_Bcast(&(grid->work_rank), 1, MPI_INT, grid->cur_rank, MPI_COMM_WORLD);
            MPI_Bcast(&(grid->gl_coord), 3, MPI_INT, grid->cur_rank, MPI_COMM_WORLD);
            /* let know the lexicographical index to other processors of next point */
            MPI_Bcast(&nl, 1, MPI_INT, (grid->cur_rank), MPI_COMM_WORLD);
            grid->cur_rank = grid->work_rank;
        }
    }
}

int main(int argc, char *argv[])
{
    /* Initialize the MPI environment */
    MPI_Init(&argc, &argv);
    /* Get the rank and size of the process */
    int my_rank, all_p;
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &all_p);

    /* Seed the random number generator */
    srand(time(NULL));

    /* Grid and array setup */
    struct GRID_INFO_T *grid = malloc(sizeof(struct GRID_INFO_T));

    Setup_grid(grid);
    define_ipt(grid);
    define_iup(grid);
    define_idn(grid);
    define_map(grid);
    /* Random walk call */
    Random_Walk(grid);

    MPI_Barrier(MPI_COMM_WORLD);
    free(grid->map);
    free(grid);
    
    /* Finalize the MPI environment */
    MPI_Finalize();
    return 0;
}