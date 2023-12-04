#include <iostream>
#include <vector>
#include <cmath>

struct SingleGrid {
    int N;
    double x0, x1, y0, y1, z0, z1;
    double ***u;
    double ***s;
    double ***res;
    double dx, dy, dz;
};

void init_grid(SingleGrid *grid, int N, double x0, double x1, double y0, double y1, double z0, double z1) {
    grid -> N = N;
    grid -> x0 = x0;
    grid -> x1 = x1;
    grid -> y0 = y0;
    grid -> y1 = y1;
    grid -> z0 = z0;
    grid -> z1 = z1;
    grid -> dx = (x1 - x0) / N;
    grid -> dy = (y1 - y0) / N;
    grid -> dz = (z1 - z0) / N;
    grid -> u = new double **[N + 1];
    grid -> s = new double **[N + 1];
    grid -> res = new double **[N + 1];
    for (int x = 0; x <= N; ++x) {
        grid -> u[x] = new double *[N + 1];
        grid -> s[x] = new double *[N + 1];
        grid -> res[x] = new double *[N + 1];
        for (int y = 0; y <= N; ++y) {
            grid -> u[x][y] = new double[N + 1];
            grid -> s[x][y] = new double[N + 1];
            grid -> res[x][y] = new double[N + 1];
        }
    }
    for (int x = 0; x <= N; ++x) {
        for (int y = 0; y <= N; ++y) {
            for (int z = 0; z <= N; ++z) {
                grid -> u[x][y][z] = 0;
                grid -> s[x][y][z] = 0;
                grid -> res[x][y][z] = 0;
            }
        }
    }
}

void del_grid(SingleGrid *grid) {
    int N = grid -> N;
    for (int x = 0; x <= N; ++x) {
        for (int y = 0; y <= N; ++y) {
            delete[] grid -> u[x][y];
            delete[] grid -> s[x][y];
            delete[] grid -> res[x][y];
        }
        delete[] grid -> u[x];
        delete[] grid -> s[x];
        delete[] grid -> res[x];
    }
    delete[] grid -> u;
    delete[] grid -> s;
    delete[] grid -> res;
}

//SingleGrid* new_grid(int N, double x0, double x1, double y0, double y1, double z0, double z1) {
//    SingleGrid *grid = new SingleGrid;
//    grid -> N = N;
//    grid -> x0 = x0;
//    grid -> x1 = x1;
//    grid -> y0 = y0;
//    grid -> y1 = y1;
//    grid -> z0 = z0;
//    grid -> z1 = z1;
//    grid -> dx = (x1 - x0) / N;
//    grid -> dy = (y1 - y0) / N;
//    grid -> dz = (z1 - z0) / N;
//    grid -> u = new double **[N + 1];
//    grid -> s = new double **[N + 1];
//    grid -> res = new double **[N + 1];
//    for (int x = 0; x <= N; ++x) {
//        grid -> u[x] = new double *[N + 1];
//        grid -> s[x] = new double *[N + 1];
//        grid -> res[x] = new double *[N + 1];
//        for (int y = 0; y <= N; ++y) {
//            grid -> u[x][y] = new double[N + 1];
//            grid -> s[x][y] = new double[N + 1];
//            grid -> res[x][y] = new double[N + 1];
//        }
//    }
//    for (int x = 0; x <= N; ++x) {
//        for (int y = 0; y <= N; ++y) {
//            for (int z = 0; z <= N; ++z) {
//                grid -> u[x][y][z] = 0;
//                grid -> s[x][y][z] = 0;
//                grid -> res[x][y][z] = 0;
//            }
//        }
//    }
//    return grid;
//}

class Multigrid {
private:
    int depth, iter;
    bool verbose;
    SingleGrid *Grids;

//    std::vector<SingleGrid> Grids;

public:
    Multigrid(int depth, int iter, double r, bool verbose);
    ~Multigrid();
    void relax(int depth);
    void relax_rb(int depth);
    void prolongation(int depth);
    void restriction(int depth);
    void residual(int depth);
    void multigrid(int lvs);
};

Multigrid::Multigrid(int depth, int iter, double r, bool verbose) {
    this -> depth = depth;
    this -> iter = iter;
    this -> verbose = verbose;
    this -> Grids = new SingleGrid[depth];
    for (int i = 0; i < depth; ++i) {
        int N = 2 << i;
        init_grid(&Grids[i], N, -r, r, -r, r, -r, r);
//        this -> Grids[i] = new_grid(N, -r, r, -r, r, -r, r);
    }
}

Multigrid::~Multigrid() {
    for (int i = 0; i < depth; ++i) {
        del_grid(&Grids[i]);
    }
    delete[] Grids;
}

void Multigrid::relax(int depth) {
    for (int i = 0; i < this -> iter; ++i) {
        this ->relax_rb(depth);
    }
}

void Multigrid::relax_rb(int depth) {
    int N = this -> Grids[depth].N;
    for (int x = 1; x < N; ++x) {
        for (int y = 1; y < N; ++y) {
            for (int z = 1 + (x + y) % 2; z < N; ++z) {
                this -> Grids[depth].u[x][y][z] = \
                (this -> Grids[depth].u[x + 1][y][z] + this -> Grids[depth].u[x - 1][y][z] +\
                this -> Grids[depth].u[x][y + 1][z] + this -> Grids[depth].u[x][y - 1][z] +\
                this -> Grids[depth].u[x][y][z + 1] + this -> Grids[depth].u[x][y][z - 1] \
                - Grids[depth].dx * Grids[depth].dx * Grids[depth].s[x][y][z]) / 6;
            }
        }
    }
    for (int x = 1; x < N; ++x) {
        for (int y = 1; y < N; ++y) {
            for (int z = 1 + (x + y + 1) % 2; z < N; ++z) {
                this -> Grids[depth].u[x][y][z] = \
                (this -> Grids[depth].u[x + 1][y][z] + this -> Grids[depth].u[x - 1][y][z] +\
                this -> Grids[depth].u[x][y + 1][z] + this -> Grids[depth].u[x][y - 1][z] +\
                this -> Grids[depth].u[x][y][z + 1] + this -> Grids[depth].u[x][y][z - 1] \
                - Grids[depth].dx * Grids[depth].dx * Grids[depth].s[x][y][z]) / 6;
            }
        }
    }

    int x = 0;
    for (int y = 1; y < N; ++y) {
        for (int z = 1; z < N; ++z) {
            int nx = 1, ny = y, nz = z;
            double r0 = sqrt((Grids[depth].x0 + Grids[depth].dx * x) *\
                             (Grids[depth].x0 + Grids[depth].dx * x) +\
                             (Grids[depth].y0 + Grids[depth].dy * y) *\
                                     (Grids[depth].y0 + Grids[depth].dy * y) +\
                             (Grids[depth].z0 + Grids[depth].dz * z) *\
                                     (Grids[depth].z0 + Grids[depth].dz * z));
            double r1 = sqrt((Grids[depth].x0 + Grids[depth].dx * nx) *\
                             (Grids[depth].x0 + Grids[depth].dx * nx) +\
                             (Grids[depth].y0 + Grids[depth].dy * ny) *\
                                     (Grids[depth].y0 + Grids[depth].dy * ny) +\
                             (Grids[depth].z0 + Grids[depth].dz * nz) *\
                                     (Grids[depth].z0 + Grids[depth].dz * nz));
            this -> Grids[depth].u[x][y][z] = Grids[depth].u[nx][ny][nz] * r1 / r0;
        }
    }

    x = N;
    for (int y = 1; y < N; ++y) {
        for (int z = 1; z < N; ++z) {
            int nx = x - 1, ny = y, nz = z;
            double r0 = sqrt((Grids[depth].x0 + Grids[depth].dx * x) *\
                             (Grids[depth].x0 + Grids[depth].dx * x) +\
                             (Grids[depth].y0 + Grids[depth].dy * y) *\
                                     (Grids[depth].y0 + Grids[depth].dy * y) +\
                             (Grids[depth].z0 + Grids[depth].dz * z) *\
                                     (Grids[depth].z0 + Grids[depth].dz * z));
            double r1 = sqrt((Grids[depth].x0 + Grids[depth].dx * nx) *\
                             (Grids[depth].x0 + Grids[depth].dx * nx) +\
                             (Grids[depth].y0 + Grids[depth].dy * ny) *\
                                     (Grids[depth].y0 + Grids[depth].dy * ny) +\
                             (Grids[depth].z0 + Grids[depth].dz * nz) *\
                                     (Grids[depth].z0 + Grids[depth].dz * nz));
            this -> Grids[depth].u[x][y][z] = Grids[depth].u[nx][ny][nz] * r1 / r0;
        }
    }

    int y = 0;
    for (int x = 0; x <= N; ++x) {
        for (int z = 1; z < N; ++z) {
            int nx = x, ny = y + 1, nz = z;
            double r0 = sqrt((Grids[depth].x0 + Grids[depth].dx * x) *\
                             (Grids[depth].x0 + Grids[depth].dx * x) +\
                             (Grids[depth].y0 + Grids[depth].dy * y) *\
                                     (Grids[depth].y0 + Grids[depth].dy * y) +\
                             (Grids[depth].z0 + Grids[depth].dz * z) *\
                                     (Grids[depth].z0 + Grids[depth].dz * z));
            double r1 = sqrt((Grids[depth].x0 + Grids[depth].dx * nx) *\
                             (Grids[depth].x0 + Grids[depth].dx * nx) +\
                             (Grids[depth].y0 + Grids[depth].dy * ny) *\
                                     (Grids[depth].y0 + Grids[depth].dy * ny) +\
                             (Grids[depth].z0 + Grids[depth].dz * nz) *\
                                     (Grids[depth].z0 + Grids[depth].dz * nz));
            this -> Grids[depth].u[x][y][z] = Grids[depth].u[nx][ny][nz] * r1 / r0;
        }
    }

    y = N;
    for (int x = 0; x <= N; ++x) {
        for (int z = 1; z < N; ++z) {
            int nx = x, ny = y - 1, nz = z;
            double r0 = sqrt((Grids[depth].x0 + Grids[depth].dx * x) *\
                             (Grids[depth].x0 + Grids[depth].dx * x) +\
                             (Grids[depth].y0 + Grids[depth].dy * y) *\
                                     (Grids[depth].y0 + Grids[depth].dy * y) +\
                             (Grids[depth].z0 + Grids[depth].dz * z) *\
                                     (Grids[depth].z0 + Grids[depth].dz * z));
            double r1 = sqrt((Grids[depth].x0 + Grids[depth].dx * nx) *\
                             (Grids[depth].x0 + Grids[depth].dx * nx) +\
                             (Grids[depth].y0 + Grids[depth].dy * ny) *\
                                     (Grids[depth].y0 + Grids[depth].dy * ny) +\
                             (Grids[depth].z0 + Grids[depth].dz * nz) *\
                                     (Grids[depth].z0 + Grids[depth].dz * nz));
            this -> Grids[depth].u[x][y][z] = Grids[depth].u[nx][ny][nz] * r1 / r0;
        }
    }

    int z = 0;
    for (int x = 0; x <= N; ++x) {
        for (int y = 0; y <= N; ++y) {
            int nx = x, ny = y, nz = z + 1;
            double r0 = sqrt((Grids[depth].x0 + Grids[depth].dx * x) *\
                             (Grids[depth].x0 + Grids[depth].dx * x) +\
                             (Grids[depth].y0 + Grids[depth].dy * y) *\
                                     (Grids[depth].y0 + Grids[depth].dy * y) +\
                             (Grids[depth].z0 + Grids[depth].dz * z) *\
                                     (Grids[depth].z0 + Grids[depth].dz * z));
            double r1 = sqrt((Grids[depth].x0 + Grids[depth].dx * nx) *\
                             (Grids[depth].x0 + Grids[depth].dx * nx) +\
                             (Grids[depth].y0 + Grids[depth].dy * ny) *\
                                     (Grids[depth].y0 + Grids[depth].dy * ny) +\
                             (Grids[depth].z0 + Grids[depth].dz * nz) *\
                                     (Grids[depth].z0 + Grids[depth].dz * nz));
            this -> Grids[depth].u[x][y][z] = Grids[depth].u[nx][ny][nz] * r1 / r0;
        }
    }

    z = N;
    for (int x = 0; x <= N; ++x) {
        for (int y = 0; y <= N; ++y) {
            int nx = x, ny = y, nz = z - 1;
            double r0 = sqrt((Grids[depth].x0 + Grids[depth].dx * x) *\
                             (Grids[depth].x0 + Grids[depth].dx * x) +\
                             (Grids[depth].y0 + Grids[depth].dy * y) *\
                                     (Grids[depth].y0 + Grids[depth].dy * y) +\
                             (Grids[depth].z0 + Grids[depth].dz * z) *\
                                     (Grids[depth].z0 + Grids[depth].dz * z));
            double r1 = sqrt((Grids[depth].x0 + Grids[depth].dx * nx) *\
                             (Grids[depth].x0 + Grids[depth].dx * nx) +\
                             (Grids[depth].y0 + Grids[depth].dy * ny) *\
                                     (Grids[depth].y0 + Grids[depth].dy * ny) +\
                             (Grids[depth].z0 + Grids[depth].dz * nz) *\
                                     (Grids[depth].z0 + Grids[depth].dz * nz));
            this -> Grids[depth].u[x][y][z] = Grids[depth].u[nx][ny][nz] * r1 / r0;
        }
    }

}

void Multigrid::prolongation(int depth) {
    int N0 = this -> Grids[depth - 1].N;
    int dx4[8] = {1, 1, 1, 1, -1, -1, -1, -1};
    int dy4[8] = {1, -1, 1, -1, 1, -1, 1, -1};
    int dz4[8] = {1, 1, -1, -1, 1, 1, -1, -1};
    for (int x = 0; x <= N0; ++x) {
        for (int y = 0; y <= N0; ++y) {
            for (int z = 0; z <= N0; ++z) {
                double u0 = Grids[depth - 1].u[x][y][z];
                // Corner
                // Edge
                // Surface

                // Inner
                if (0 < x && x < N0 && 0 < y && y < N0 && 0 < z && z < N0) {
                    for (int ind = 0; ind < 8; ++ind) {
                        this -> Grids[depth].u[(x << 1) + dx4[ind]][(y << 1) + dy4[ind]][(z << 1) + dz4[ind]] += u0 / 8;
                    }
                }
                else if ()
            }
        }
    }

//    for (int x = 0; x <= this -> Grids[depth].N; ++x) {
//        for (int y = 0; y <= this -> Grids[depth].N; ++y) {
//            for (int z = 0; z <= this -> Grids[depth].N; ++z) {
//                if ((x & 1) || (y & 1) || (z & 1)) {
//
//                }
//                else {
//                    this -> Grids[depth].u[x][y][z] = this -> Grids[depth - 1].u[x >> 1][y >> 2][z >> 2];
//                }
//            }
//        }
//    }
}

void Multigrid::restriction(int depth) {

}

void Multigrid::residual(int depth) {

}

void Multigrid::multigrid(int lvs) {

}

int main() {
    Multigrid mgs = Multigrid(5, 1, 10, false);

    return 0;
}
