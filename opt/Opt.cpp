#include <iostream>
#include <math.h>
#include <random>
#include <fstream>

/* basic physical constants */
#define SIGMA 1
#define EP 1
#define M 1
#define R_C 2.5
#define DELTA 0.01631689
#define MAX_DISPLACE 0.1

void init_coord(float **coord, int num_particle, float size_box) {
    /* generate random number betwwen 0 and 1 */
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<float> dist_r(0.0,1.0);
    /* initiate coordination matrix with random numbers */
//    float gap = size_box/float(pow(num_particle,0.33333));
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            coord[i][j] = dist_r(gen) * size_box;
        }
    }
}


void init_0_matrix(float **matrix, int num_particle) {
    /* initiate matrix with 0 */
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            matrix[i][j] = 0.0;
        }
    }
}


float sum_matrix(float **matrix, int num_particle) {
    /* sum up all elements */
    float sum = 0.0;
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            sum += fabs(matrix[i][j]);
        }
    }
    return sum;
}


void print_coord(float **coord, float num_particle) {
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            std::cout << coord[i][j] << "\t";
        }
        std::cout << "\n";
    }
    std::cout << "\n";
}


float min_image(float coord_i, float coord_j, float size_box) {
    float a = coord_j-size_box-coord_i;
    float b = coord_j-coord_i;
    float c = coord_j+size_box-coord_i;
    float abs_a = fabs(a);
    float abs_b = fabs(b);
    float abs_c = fabs(c);
    if (abs_a < abs_b) {
        if (abs_a < abs_c) {
            return a;
        }
        else {
            return c;
        }
    }
    else {
        if (abs_b < abs_c) {
            return b;
        }
        else {
            return c;
        }
    }
}


void calcu_force(float **coord, float **force, int num_particle, float size_box, float &V) {
    float r, k, k2, F;
    float iract[num_particle][num_particle][3];
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < num_particle; j++) {
            iract[i][j][0] = 0.0;
            iract[i][j][1] = 0.0;
            iract[i][j][2] = 0.0;
        }
    }

    for (int i = 0; i < num_particle; i++) {
        for (int j = i + 1; j < num_particle; j++) {
            float dx, dy, dz;
            /* check of nearest mirror image */
            dx = min_image(coord[i][0], coord[j][0], size_box);
            dy = min_image(coord[i][1], coord[j][1], size_box);
            dz = min_image(coord[i][2], coord[j][2], size_box);
            /* check if interaction exists */
            r = sqrt(dx * dx + dy * dy + dz * dz);
//            std::cout << dx << "\t" << dy << "\t" << dz << "\t" << r << std::endl;
            if (r < R_C) {
                k = pow((SIGMA/r),6);
                k2 = k * k;
                /* potential energy */
                V += 4 * EP * (k2 - k) + DELTA;
                /* force */
                F = 24 * EP * (-2 * k2/r + k/r);
 //               std::cout << F << std::endl;
                /* add forces by other atoms within R_C of three directions */
                /* Newton third law */
                iract[i][j][0] = F * dx/r;
                iract[j][i][0] = -iract[i][j][0];
                iract[i][j][1] = F * dy/r;
                iract[j][i][1] = -iract[i][j][1];
                iract[i][j][2] = F * dz/r;
                iract[j][i][2] = -iract[i][j][2];
            }
        }
    }
    /* add all forces together for three directions */
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < num_particle; j++) {
            force[i][0] += iract[i][j][0];
            force[i][1] += iract[i][j][1];
            force[i][2] += iract[i][j][2];
        }
    }
}


float boundary_coord(float coord, float size_box) {
    if (coord < 0) {
        return coord + size_box;
    }
    else if (coord > size_box) {
        return coord - size_box;
    }
    else {
        return coord;
    }
}


void calcu_all(float **coord, float **veloc, float **force, float **accle, int num_particle, float size_box, float dt, float &V) {
    
    /* initiate acceleration matrix */
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            accle[i][j] = force[i][j]/M;
        }
    }
    
    /* calculate new coordination matrix */
    float dcoord;
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            dcoord = veloc[i][j] * dt + 0.5 * accle[i][j] * accle[i][j];
            if (fabs(dcoord) > MAX_DISPLACE) {
                dcoord = ((dcoord>0)?1:(-1)) * MAX_DISPLACE;
            }
            coord[i][j] = boundary_coord(coord[i][j] + dcoord, size_box);
        }
    }
    
    /* calculate new velocity matrix step 1 */
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            veloc[i][j] = 0.5 * dt * accle[i][j];
        }
    }
    
    /* calculate new force matrix */
    calcu_force(coord, force, num_particle, size_box, V);
    
    /* transfer force matrix to acceleration matrix */
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            accle[i][j] = force[i][j]/M;
        }
    }
    
    /* calculate new velocity matrix step 2 */
    for (int i = 0; i < num_particle; i++) {
        for (int j = 0; j < 3; j++) {
            veloc[i][j] += 0.5 * dt * accle[i][j];
        }
    }
}


int main() {
    
    int num_particle;
    std::cout << "Input number of particle:" << std::endl;
    std::cin >> num_particle;
    
    float size_box;
    std::cout << "Input size of box:" << std::endl;
    std::cin >> size_box;
    
    int num_md;
    std::cout << "Input maximum iteration number of MD:" << std::endl;
    std::cin >> num_md;
    
    float dt = 0.2 * size_box/float(pow(num_particle, 0.33333));
    float dt_save = dt;
    
    
    /* coordination matrix of num_particle */
    float **coord;
    coord = (float**)malloc(num_particle*sizeof(float*));
    for (int n = 0; n < num_particle; n++) {
        coord[n] =(float*)malloc(3*sizeof(float));
    }
    
    /* force matrix of num_particle */
    float **force;
    force = (float**)malloc(num_particle*sizeof(float*));
    for (int n = 0; n < num_particle; n++) {
        force[n] =(float*)malloc(3*sizeof(float));
    }
    
    /* acceleration matrix of num_particle */
    float **accle;
    accle = (float**)malloc(num_particle*sizeof(float*));
    for (int n = 0; n < num_particle; n++) {
        accle[n] =(float*)malloc(3*sizeof(float));
    }
    
    /* velocity matrix of num_particle */
    float **veloc;
    veloc = (float**)malloc(num_particle*sizeof(float*));
    for (int n = 0; n < num_particle; n++) {
        veloc[n] =(float*)malloc(3*sizeof(float));
    }
    
    /* initiate coordination and velocity matrix */
    init_coord(coord, num_particle, size_box);
    init_0_matrix(veloc, num_particle);
    
    /* xyz output prepare */
    std::ofstream xyz_file;
    xyz_file.open("md.xyz");
    /* std output prepare */
    std::cout << "iterations" << "\t"  << "   " << "force" << "\t" << "   " << "potential energy" << std::endl;
    
    /* main calcu */
    float V, V_old = -65535.0, F[num_md];
    for (int iter = 0; iter < num_md; iter++) {
        V = 0;
        /* initiate force matrix */
        init_0_matrix(force, num_particle);
        
        /* calculate new coordination matrix and velocity matrix */
        calcu_all(coord, veloc, force, accle, num_particle, size_box, dt, V);
        
        /* print info of this md round */
        /* xyz output */
        xyz_file << num_particle << "\n" << "t = " << iter * dt << std::endl;
        for (int i = 0; i < num_particle; i++) {
            xyz_file << "   " << "Ar" << "   " << coord[i][0] << "   " << coord[i][1] << "   " << coord[i][2] << std::endl;
        }
        
        /* std output */
        F[iter] = sum_matrix(force, num_particle);
        std::cout << "    " << "No." << iter << "\t" << F[iter] << "\t" << 2 * V << "\t" << dt << std::endl;
        
        /* check force and potential energy differences */
        if ((F[iter]/num_particle < 0.1) && ((V - V_old)/num_particle < 0.0001)) {
            break;
        }
        
        /* aviod vibration between several point */
        if (F[iter] > 1.1 * F[iter-1]) {
            dt *= 0.9;
        }
        /* aviod slow convergence */
        else if (dt < dt_save && F[iter] > 0.9 * F[iter-1]) {
            dt *= 1.1;
        }
        
        V_old = V;
    }
    
    xyz_file.close();
}
