#ifndef __RANDWERTE_H__
#define __RANDWERTE_H__


/**
 * The boundary values of the problem are set. Top Wall Velocity is 1.

 Note: Assume for example: imax = 4, jmax = 4, the standard to be followed is:
    ->  The domain for each variable is constructed ALONGWITH the ghost cells.
        Thus, the entry for the matrix() and free_matrix() functions would be:

            matrix(0, imax+1, 0, jmax+1);

        and

            free_matrix(0, imax+1, 0, jmax+1);

    ->  The matrix is to be organised as ROW-MAJOR, i.e. the data starts from
        (i=0,j=0) and proceeds as (i=1,j=0) (i=2,j=0) and so on.

    ->  In the following diagram, the data is stored in the CELL-CENTERS, and
        NOT on the nodes. :)

    +------+------+------+------+------+------+
    |      |      |      |      |      |      |
    |      |      |      |      |      |      |
    +------+------+------+------+------+------+
    |      |      |      |      |      |      |
    |      |      |      |      |      |      |
    +------+------+------+------+------+------+
    |      |      |      |      |      |      |
    |      |      |      |      |      |      |
    +------+------+------+------+------+------+
    |  .   |      |   .  |      |      |      |
    |  .   |      | .    |      |      |      |
    +------+------+------+------+------+------+
    |      |      |      |      |      |      |
    | 0,1  |      |      |      |      |      |
    +------+------+------+------+------+------+
    |      |      |      |      |      |      |
    | 0,0  | 1,0  | ..   |      |      |      |
    +------+------+------+------+------+------+


 */
void boundaryvalues(
  int imax,
  int jmax,
  double **U,
  double **V
);

#endif
