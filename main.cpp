#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include <fstream>
#include <tuple>

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


class Multigrid {
private:
    int depth, iter, maxIter;
    bool verbose;
    SingleGrid *Grids;
    double ***invalpha;
    double ***alpha;
    double ***beta;
    double ***init_u;
    double TP_epsilon = 1E-6;
    double TP_tiny = 0;
//    double ***res;
//    std::vector<SingleGrid> Grids;

public:
    Multigrid(int depth, int iter, double r, bool verbose, int maxIter);
    ~Multigrid();
    void relax(int depth);
    void relax_rb(int depth);
    void prolongation(int depth);
    void restriction(int depth);
    void residual(int depth);
    void multigrid(int lvs);
    double norm_residual();
    void solve_lin();
    void solve();
    void init(const std::string& path_name);
    void write_u(const std::string& path_name);
    void write_psi(const std::string& path_name);
};

Multigrid::Multigrid(int depth, int iter, double r, bool verbose, int maxIter) {
    this -> depth = depth;
    this -> iter = iter;
    this -> verbose = verbose;
    this -> Grids = new SingleGrid[depth];
    this -> maxIter = maxIter;
    for (int i = 0; i < depth; ++i) {
        int N = 2 << i;
        init_grid(&Grids[i], N, -r, r, -r, r, -r, r);
//        this -> Grids[i] = new_grid(N, -r, r, -r, r, -r, r);
    }
    int N = 1 << depth;
    this -> invalpha = new double **[N + 1];
    this -> alpha = new double **[N + 1];
    this -> beta = new double **[N + 1];
    this -> init_u = new double **[N + 1];
//    this -> res = new double **[N + 1];
    for (int x = 0; x <= N; ++x) {
        this -> invalpha[x] = new double *[N + 1];
        this -> alpha[x] = new double *[N + 1];
        this -> beta[x] = new double *[N + 1];
        this -> init_u[x] = new double *[N + 1];
//        this -> res[x] = new double *[N + 1];
        for (int y = 0; y <= N; ++y) {
            this -> invalpha[x][y] = new double [N + 1];
            this -> alpha[x][y] = new double [N + 1];
            this -> beta[x][y] = new double [N + 1];
            this -> init_u[x][y] = new double [N + 1];
//            this -> res[x][y] = new double [N + 1];
            for (int z = 0; z <= N; ++z) {
                this -> init_u[x][y][z] = 0;
            }
        }
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
                - this -> Grids[depth].dx * this -> Grids[depth].dx * this -> Grids[depth].s[x][y][z]) / 6;
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
                - this -> Grids[depth].dx * this -> Grids[depth].dx * this -> Grids[depth].s[x][y][z]) / 6;
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
    int dx8[8] = {1, 1, 1, 1, -1, -1, -1, -1};
    int dy8[8] = {1, -1, 1, -1, 1, -1, 1, -1};
    int dz8[8] = {1, 1, -1, -1, 1, 1, -1, -1};
    int dx4[4] = {1, 1, -1, -1};
    int dy4[4] = {1, -1, 1, -1};
    int dx2[2] = {1, -1};
    for (int x = 0; x <= N0; ++x) {
        for (int y = 0; y <= N0; ++y) {
            for (int z = 0; z <= N0; ++z) {
                double u0 = Grids[depth - 1].u[x][y][z];
                // Inner
                if (0 < x && x < N0 && 0 < y && y < N0 && 0 < z && z < N0) {
                    for (int ind = 0; ind < 8; ++ind) {
                        this -> Grids[depth].u[(x << 1) + dx8[ind]][(y << 1) + dy8[ind]][(z << 1) + dz8[ind]] += u0 / 8;
                    }
                }
                // Surface
                else if (0 < x && x < N0 && 0 < y && y < N0) {
                    for (int ind = 0; ind < 4; ++ind) {
                        this -> Grids[depth].u[(x << 1) + dx4[ind]][(y << 1) + dy4[ind]][(z << 1)] += u0 / 4;
                    }
                }
                else if (0 < x && x < N0 && 0 < z && z < N0) {
                    for (int ind = 0; ind < 4; ++ind) {
                        this -> Grids[depth].u[(x << 1) + dx4[ind]][(y << 1)][(z << 1) + dy4[ind]] += u0 / 4;
                    }
                }
                else if (0 < z && z < N0 && 0 < y && y < N0) {
                    for (int ind = 0; ind < 4; ++ind) {
                        this -> Grids[depth].u[(x << 1)][(y << 1) + dy4[ind]][(z << 1)  + dx4[ind]] += u0 / 4;
                    }
                }
                // Edge
                else if (0 < x && x < N0) {
                    for (int ind = 0; ind < 2; ++ind) {
                        this -> Grids[depth].u[(x << 1) + dx2[ind]][(y << 1)][(z << 1)] += u0 / 2;
                    }
                }
                else if (0 < y && y < N0) {
                    for (int ind = 0; ind < 2; ++ind) {
                        this -> Grids[depth].u[(x << 1)][(y << 1) + dx2[ind]][(z << 1)] += u0 / 2;
                    }
                }
                else if (0 < z && z < N0) {
                    for (int ind = 0; ind < 2; ++ind) {
                        this -> Grids[depth].u[(x << 1)][(y << 1)][(z << 1) + dx2[ind]] += u0 / 2;
                    }
                }
                //Corner
                else {
                    this -> Grids[depth].u[x << 1][y << 1][z << 1] = u0;
                }
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
    int N0 = this -> Grids[depth - 1].N;
    int dx16[6] = {0, 0, 0, 0, -1, 1};
    int dy16[6] = {0, 1, 0, -1, 0, 0};
    int dz16[6] = {1, 0, -1, 0, 0, 0};
    int dx64[8] = {1, 1, 1, 1, -1, -1, -1, -1};
    int dy64[8] = {1, -1, 1, -1, 1, -1, 1, -1};
    int dz64[8] = {1, 1, -1, -1, 1, 1, -1, -1};
    int dx32[12] = {0, 1, 0, -1, 1, 1, -1, -1, 0, 1, 0, -1};
    int dy32[12] = {1, 0, -1, 0, 1, -1, -1, 1, 1, 0, -1, 0};
    int dz32[12] = {1, 1, 1, 1, 0, 0, 0, 0, -1, -1, -1, -1};

    int dx2_16[4] = {1, 1, -1, -1};
    int dy2_16[4] = {1, -1, 1, -1};
    int dx2_8[4] = {1, 0, -1, 0};
    int dy2_8[4] = {0, -1, 0, 1};
    for (int x = 0; x <= N0; ++x) {
        for (int y = 0; y <= N0; ++y) {
            for (int z = 0; z <= N0; ++z) {
                // Inner
                if (0 < x && x < N0 && 0 < y && y < N0 && 0 < z && z < N0) {
                    this -> Grids[depth].s[x][y][z] = this -> Grids[depth + 1].res[x << 1][y << 1][z << 1] / 8;
                    for (int ind = 0; ind < 6; ++ind) {
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1) + dx16[ind]][(y << 1) + dy16[ind]][(z << 1) + dz16[ind]] / 16;
                    }
                    for (int ind = 0; ind < 12; ++ind) {
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1) + dx32[ind]][(y << 1) + dy32[ind]][(z << 1) + dz32[ind]] / 32;
                    }
                    for (int ind = 0; ind < 8; ++ind) {
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1) + dx64[ind]][(y << 1) + dy64[ind]][(z << 1) + dz64[ind]] / 64;
                    }
                }

                // Surface
                else if (0 < x && x < N0 && 0 < y && y < N0) {
                    this -> Grids[depth].s[x][y][z] = this -> Grids[depth + 1].res[x << 1][y << 1][z << 1] / 4;
                    for (int ind = 0; ind < 4; ++ind) {
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1) + dx2_16[ind]][(y << 1) + dy2_16[ind]][(z << 1)] / 16;
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1) + dx2_8[ind]][(y << 1) + dy2_8[ind]][(z << 1)] / 8;
                    }
                }
                else if (0 < x && x < N0 && 0 < z && z < N0) {
                    this -> Grids[depth].s[x][y][z] = this -> Grids[depth + 1].res[x << 1][y << 1][z << 1] / 4;
                    for (int ind = 0; ind < 4; ++ind) {
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1) + dx2_16[ind]][(y << 1)][(z << 1) + dy2_16[ind]] / 16;
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1) + dx2_8[ind]][(y << 1)][(z << 1) + dy2_8[ind]] / 8;
                    }
                }
                else if (0 < z && z < N0 && 0 < y && y < N0) {
                    this -> Grids[depth].s[x][y][z] = this -> Grids[depth + 1].res[x << 1][y << 1][z << 1] / 4;
                    for (int ind = 0; ind < 4; ++ind) {
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1)][(y << 1) + dx2_16[ind]][(z << 1) + dy2_16[ind]] / 16;
                        this -> Grids[depth].s[x][y][z] += this -> Grids[depth + 1].res[(x << 1)][(y << 1) + dx2_8[ind]][(z << 1) + dy2_8[ind]] / 8;
                    }
                }

                // Edge
                else if (0 < x && x < N0) {
                    this -> Grids[depth].s[x][y][z] = this -> Grids[depth + 1].res[x << 1][y << 1][z << 1] / 2 + this -> Grids[depth + 1].res[(x << 1) - 1][y << 1][z << 1] / 4 + this -> Grids[depth + 1].res[(x << 1) + 1][y << 1][z << 1] / 4;
                }
                else if (0 < y && y < N0) {
                    this -> Grids[depth].s[x][y][z] = this -> Grids[depth + 1].res[x << 1][y << 1][z << 1] / 2 + this -> Grids[depth + 1].res[(x << 1)][(y << 1) - 1][z << 1] / 4 + this -> Grids[depth + 1].res[(x << 1)][(y << 1) + 1][z << 1] / 4;
                }
                else if (0 < z && z < N0) {
                    this -> Grids[depth].s[x][y][z] = this -> Grids[depth + 1].res[x << 1][y << 1][z << 1] / 2 + this -> Grids[depth + 1].res[(x << 1)][y << 1][(z << 1) - 1] / 4 + this -> Grids[depth + 1].res[(x << 1)][y << 1][(z << 1) + 1] / 4;
                }
                //Corner
                else {
                    this -> Grids[depth].s[x][y][z] = this -> Grids[depth + 1].res[x << 1][y << 1][z << 1];
                }
            }
        }
    }

//    for (int x = 0; x <= N0; ++x) {
//        for (int y = 0; y <= N0; ++y) {
//            for (int z = 0; z <= N0; ++z) {
//                this -> Grids[depth].u[x][y][z] = 0;
//            }
//        }
//    }
}

void Multigrid::residual(int depth) {
    int N0 = this -> Grids[depth - 1].N;
    for (int x = 1; x < N0; ++x) {
        for (int y = 1; y < N0; ++y) {
            for (int z = 1; z < N0; ++z) {
                this -> Grids[depth].res[x][y][z] = this -> Grids[depth].s[x][y][z] - (this -> Grids[depth].u[x + 1][y][z] + this -> Grids[depth].u[x - 1][y][z] +\
                this -> Grids[depth].u[x][y + 1][z] + this -> Grids[depth].u[x][y - 1][z] +\
                this -> Grids[depth].u[x][y][z + 1] + this -> Grids[depth].u[x][y][z - 1] \
                - 6 * this -> Grids[depth].u[x][y][z]) / this -> Grids[depth].dx / this -> Grids[depth].dx;
            }
        }
    }
}

void Multigrid::multigrid(int lvs) {
    if (lvs == 1) {
        this -> relax(lvs);
        return;
    }
    this -> relax(lvs);
    this -> residual(lvs);
    this -> restriction(lvs - 1);
    this -> multigrid(lvs - 1);
    this -> prolongation(lvs);
    this -> relax(lvs);
}

void Multigrid::solve_lin() {
    if (this -> verbose) {
        std::cout << "Solving Linear PDE...\n";
    }
    while (true) {
        this -> multigrid(this -> depth - 1);
        double x = this -> norm_residual();
        if (this -> verbose) {
//            std::cout << "  local norm residual: " << x << '\r';
//            printf("  local norm residual: %f\r", x);
        }
        if (x < 1E-12) {
            std::cout << "\n";
            break;
        }
    }
    if (this -> verbose) {
        std::cout << "Done!\n";
    }
}

double Multigrid::norm_residual() {
    double cnt = 0;
    int N0 = this -> Grids[this -> depth - 1].N;
    double dx = this -> Grids[this -> depth - 1].dx;
    for (int x = 1; x < N0; ++x) {
        for (int y = 1; y < N0; ++y) {
            for (int z = 1; z < N0; ++z) {
                cnt += this -> Grids[this -> depth - 1].res[x][y][z] * this -> Grids[this -> depth - 1].res[x][y][z] / (N0 - 1);
            }
        }
    }

    return sqrt(cnt) * dx * dx;
}

void Multigrid::solve() {
    int N0 = this -> Grids[depth - 1].N;
    int depth = this -> depth - 1;
    for (int iter = 0; iter < this -> maxIter; ++iter) {
        double mean_s = 0;

//        for (int d = 0; d <= depth; ++d) {
//            int n = this -> Grids[d].N;
//            for (int x = 0; x <= n; ++x) {
//                for (int y = 0; y <= n; ++y) {
//                    for (int z = 0; z <= n; ++z) {
//                        this -> Grids[d].u[x][y][z] = 0;
//                    }
//                }
//            }
//        }



        for (int x = 1; x < N0; ++x) {
            for (int y = 1; y < N0; ++y) {
                for (int z = 1; z < N0; ++z) {
                    this -> Grids[depth].s[x][y][z] = - this-> beta[x][y][z] * std::pow(this -> alpha[x][y][z] * (1 + this -> init_u[x][y][z]) + 1, -7) - (this -> init_u[x + 1][y][z] + this -> init_u[x - 1][y][z] +\
                                                        this -> init_u[x][y + 1][z] + this -> init_u[x][y - 1][z] +\
                                                        this -> init_u[x][y][z + 1] + this -> init_u[x][y][z - 1] \
                                                        - 6 * this -> init_u[x][y][z]) / this -> Grids[depth].dx / this -> Grids[depth].dx;
                    mean_s += this -> Grids[depth].s[x][y][z] * this -> Grids[depth].s[x][y][z];
                }
            }
        }
        std::cout << "Mean source: " << std::sqrt(mean_s / (N0 - 1) / (N0 - 1) / (N0 - 1)) << "\n";
        double nres = 0;
        this -> solve_lin();
        for (int x = 0; x <= N0; ++x) {
            for (int y = 0; y <= N0; ++y) {
                for (int z = 0; z <= N0; ++z) {
                    this -> init_u[x][y][z] += this -> Grids[depth].u[x][y][z];
                    nres += this -> Grids[depth].u[x][y][z] * this -> Grids[depth].u[x][y][z];
                }
            }
        }
        nres = std::sqrt(nres / (N0 + 1) / (N0 + 1) / (N0 + 1));
        if (this -> verbose) {
            std::cout << "Iterating.. residual: " << nres << "\n";
        }
    }
}

void Multigrid::init(const std::string& path_name) {
    std::ifstream fin;
    fin.open(path_name);
    if (!fin) {
        std::cout << "Fail to open file.\n";
        return;
    }
    std::cout << "Setting Initial alpha, beta...\n";
    int NofP;
    fin >> NofP;
    std::vector<std::tuple<double, double, double, double, double, double, double, double, double, double>> punctures;

    for (int i = 0; i < NofP; ++i) {
        double Mi, xi, yi, zi, Pix, Piy, Piz, Six, Siy, Siz;
        fin >> Mi >> xi >> yi >> zi >> Pix >> Piy >> Piz >> Six >> Siy >> Siz;
        punctures.emplace_back(Mi, xi, yi, zi, Pix, Piy, Piz, Six, Siy, Siz);
        std::cout << "xyz: " << xi << " " << yi << " " << zi << "\n";
    }
    int depth = this -> depth - 1;
    int N0 = this -> Grids[depth].N;
    double alpha_inv, si, si2, si3, alpha_p, beta_p, Aij;
    double nx, ny, nz;
    double ix, iy, iz;
    double P[3];
    double S[3];
    double n[3];
    double n_P, n_S[3];
    for (int x = 0; x <= N0; ++x) {
        for (int y = 0; y <= N0; ++y) {
            for (int z = 0; z <= N0; ++z) {
                alpha_inv = 0;
                nx = this -> Grids[depth].x0 + this -> Grids[depth].dx * x;
                ny = this -> Grids[depth].y0 + this -> Grids[depth].dx * y;
                nz = this -> Grids[depth].z0 + this -> Grids[depth].dx * z;
//                std::cout << "nxyz: " << nx << " " << ny << " " << nz << "\n";
//                si = sqrt((this -> Grids[depth].x0 + this -> Grids[depth].dx * x) *\
//                             (this -> Grids[depth].x0 + this -> Grids[depth].dx * x) +\
//                             (this -> Grids[depth].y0 + this -> Grids[depth].dy * y) *\
//                                     (this -> Grids[depth].y0 + this -> Grids[depth].dy * y) +\
//                             (this -> Grids[depth].z0 + this -> Grids[depth].dz * z) *\
//                                     (this -> Grids[depth].z0 + this -> Grids[depth].dz * z));
//                si2 = si * si;
//                si4 = si2 * si2;
//                si = std::max(this -> TP_tiny, std::pow(si4 + std::pow(this -> TP_epsilon, 4), 1/4));
//                si2 = si * si;
//                si4 = si2 * si2;
                beta_p = 0;
                for (int ind = 0; ind < NofP; ++ind) {
                    ix = std::get<1>(punctures[ind]);
                    iy = std::get<2>(punctures[ind]);
                    iz = std::get<3>(punctures[ind]);
                    si2 = (nx - ix) * (nx - ix) + (ny - iy) * (ny - iy) + (nz - iz) * (nz - iz);
                    si2 = std::sqrt(std::pow(si2, 2) + std::pow(TP_epsilon, 4));
                    if (si2 < pow(TP_tiny, 2)) si2 = pow(TP_tiny, 2);
                    si = std::sqrt(si2);
//                    std::cout << "si: " << si << "\n";
                    alpha_inv += std::get<0>(punctures[ind]) / 2 / si;
//                    if (x == 19 && y == 16 && z == 16) {
//                        std::cout << "xyzs: " << alpha_inv << " " << si << " " << ix << ", " << iy << ", " << iz << ", " << nx << ", " << ny << ", " << nz << "\n";
//                    }
                }

                for (int i = 0; i < 3; ++i) {
                    for (int j = 0; j < 3; ++j) {
                        Aij = 0;
                        for (int ind = 0; ind < NofP; ++ind) {
                            ix = std::get<1>(punctures[ind]);
                            iy = std::get<2>(punctures[ind]);
                            iz = std::get<3>(punctures[ind]);
                            si2 = (nx - ix) * (nx - ix) + (ny - iy) * (ny - iy) + (nz - iz) * (nz - iz);
                            si2 = std::sqrt(std::pow(si2, 2) + std::pow(TP_epsilon, 4));
                            if (si2 < pow(TP_tiny, 2)) si2 = pow(TP_tiny, 2);
                            si = std::sqrt(si2);
                            si3 = si * si2;

                            P[0] = std::get<4>(punctures[ind]);
                            P[1] = std::get<5>(punctures[ind]);
                            P[2] = std::get<6>(punctures[ind]);
                            S[0] = std::get<7>(punctures[ind]);
                            S[1] = std::get<8>(punctures[ind]);
                            S[2] = std::get<9>(punctures[ind]);
                            n[0] = (nx - ix) / si;
                            n[1] = (ny - iy) / si;
                            n[2] = (nz - iz) / si;
                            n_P = 0;
                            n_S[0] = n[1] * S[2] - n[2] * S[1];
                            n_S[1] = n[2] * S[0] - n[0] * S[2];
                            n_S[2] = n[0] * S[1] - n[1] * S[0];
                            for (int k = 0; k < 3; ++k) {
                                n_P += n[k] * P[k];
                            }
                            Aij += 1.5 * (P[i] * n[j] + P[j] * n[i] + n_P * n[i] * n[j]) / si2 - 3.0 * (n_S[i] * n[j] + n_S[j] * n[i]) / si3;
                            if (i == j) {
                                Aij -= 1.5 * (n_P / si2);
                            }
                        }
                        beta_p += Aij * Aij;
                    }
                }
                alpha_p = 1 / alpha_inv;
                beta_p = beta_p * std::pow(alpha_p, 7) / 8;
                this -> invalpha[x][y][z] = alpha_inv;
                this -> alpha[x][y][z] = alpha_p;
                this -> beta[x][y][z] = beta_p;
                std::cout << nx << " " << ny << " " << nz << " " << alpha_p << " " << beta_p << "\n";
            }
        }
    }
    fin.close();
    std::cout << "  Done!\n";
}

void Multigrid::write_u(const std::string& path_name) {
    std::ofstream fout;
    fout.open(path_name);
    fout.precision(20);
    int depth = this -> depth - 1;
    int N0 = this -> Grids[depth].N;
    fout << N0 << "\n";
    for (int x = 0; x <= N0; ++x) {
        for (int y = 0; y <= N0; ++y) {
            for (int z = 0; z <= N0; ++z) {
                fout << this -> init_u[x][y][z] << " ";
            }
            fout << "\n";
        }
        fout << "\n";
    }
    fout.close();
}

void Multigrid::write_psi(const std::string &path_name) {
    std::ofstream fout;
    fout.open(path_name);
    fout.precision(20);
    int depth = this -> depth - 1;
    int N0 = this -> Grids[depth].N;
    fout << N0 << "\n";
    for (int x = 0; x <= N0; ++x) {
        for (int y = 0; y <= N0; ++y) {
            for (int z = 0; z <= N0; ++z) {
                fout << 1 + this -> invalpha[x][y][z] + this -> init_u[x][y][z] << " ";
            }
            fout << "\n";
        }
        fout << "\n";
    }
    fout.close();
}

int main() {
    Multigrid mgs = Multigrid(5, 1, 6, true, 10);
    mgs.init("data1.in");
//    mgs.solve();
//    mgs.write_u("out.txt");
//    mgs.write_psi("data-init-4-psi.out");
    return 0;
}
