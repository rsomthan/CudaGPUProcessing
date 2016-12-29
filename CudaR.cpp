#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

typedef struct { double x, y; int group; } point_t, *point;
// x and y are coordinates
// group contains the cluster number it is associated with
double randf(double m)
{
	return m * rand() / (RAND_MAX - 1.); // RAND_MAX is the max value rand() can return it is always grater than 32767
}
 // count is the number of data-points we want
point gen_xy(int count, double radius)//generating the data-points in space
{
	double ang, r;
	point p, pt = malloc(sizeof(point_t) * count);//if suppose we want 100000 data-points, so it will allocate 100000 different structures for these 100000 data-points. These arrays
	
// are pointed by the pointer pt
//at present there is nothing in p

	/* note: this is not a uniform 2-d distribution */
	for (p = pt + count; p-- > pt;) {//iterating it 'count' no. of times ; till p>pt, the loop will run
	
	
		ang = randf(2 * 3.14); //M_PI is constant whose value is pi & here ang can never be equal to 2pi
		r = randf(radius);
		p->x = r * cos(ang);
		p->y = r * sin(ang);
	}
 //after this loop p will be equal to pt(decrementing after one iteration; loop will stop at p=pt)
	return pt;//at this juncture, p is equal to pt. So returning pt and p is one and the same thing
}
 
//inline
	double dist2(point a, point b) // used to calculate euclidian distance 
{
	double x = a->x - b->x, y = a->y - b->y;
	return (x*x + y*y);
}
 
//inline 
	int nearest(point pt, point cent, int n_cluster, double *d2)
{
	int i, min_i;
	point c;
	double d, min_d;
 
#	define for_n for (c = cent, i = 0; i < n_cluster; i++, c++) //if we are currently in 5th centroid then loop will run from 0 to 4 
	for_n { // this outer loop is of no use
		min_d = HUGE_VAL; // returns a value with is not in range of double datatype 
		min_i = pt->group;
		for_n {
			if (min_d > (d = dist2(c, pt))) // dist2 calculates the distance of a data-point with all the clusters 
			{
				min_d = d; min_i = i; //i is current cluster(i.e. point nearest to this cluster) 
			}
		}
	}
	if (d2) //if d2!=0.... this is used later in do while loop
	*d2 = min_d; 
	return min_i;
}
 
void kpp(point pts, int len, point cent, int n_cent)
{
#	define for_len for (j = 0, p = pts; j < len; j++, p++) //j is going from 0 to 1000000 and p is from the addr of first data-point to the addr of last data-point
	int i, j;
	int n_cluster;
	double sum;
    double *d = malloc(sizeof(double) * len); //allocation of memory to d so that it can point to array of all data points together
	
	point p, c;
	cent[0] = pts[ rand() % len ]; //cent[0] is pointing to the first centroid, pts = v= pointer to the entire dataset
	//cent[0] is storing any one random point from the dataset as the first centroid
	for (n_cluster = 1; n_cluster <= n_cent; n_cluster++) //here n_cluster is a local variable ranging from 1 to no. of clusters
	{
		sum = 0;
		for_len {
			p->group=nearest(p, cent, n_cluster, d + j);//d+j=d[j]
			// we have changed this ^^^ line
			sum += d[j];//it is the sum of least distances of all data points with its respective centroid eg. sum=1000
		}
		sum = randf(sum); //here we are generating a random no. which is less than the sum i.e. 1000
		//bcaz 1000*random number/32767 will give a number less than 1000
		for_len {
			if ((sum -= d[j]) > 0) //here d[j] is the distance of first data point from its centroid
			continue;
			cent[n_cluster] = pts[j]; //see the dry run in himani's book
			break;
		}
	}
	//after this loop we have all centroids with a somewhat random value
	for_len p->group = nearest(p, cent, n_cluster, 0);  // <-- we have commented this bcaz not needed
	free(d);
}
 
point lloyd(point pts, int len, int n_cluster)
{
	int i, j, min_i;
	int changed;
 
	point cent = malloc(sizeof(point_t) * n_cluster), p, c;// allocating size of all centroids and pointed by cent
 
	/* assign init grouping randomly */
	//for_len p->group = j % n_cluster;
 
	/* or call k++ init */
	kpp(pts, len, cent, n_cluster);
 
	do {
		/* group element for centroids are used as counters */
		for_n 
		//for (c = cent, i = 0; i < n_cluster; i++, c++)
		{ c->group = 0; c->x = c->y = 0; } //basically clearing the structure c and not cent
		
		
		
		//problem strts
		
		for_len
		//for (j = 0, p = pts; j < len; j++, p++)
		{
		
			c = cent + p->group; //addr stored in c = addr of first centroid + group no. of first data-point i.e. addr centroid of first data point
			//read himani's note book
			//basically we are coping the group from p to c
			
			++c-> group; //bcaz we have done c->group=0 this makes cluster size from 0 to 10 but we want it to be 1 to 11 thus we are incrementing by 1
			//we have changed post to pre
			c->x += p->x; c->y += p->y;// this is the sum of all x and y coordinates of p
		}
		//after this loop only c[0] the sum of x and y coordinates of 100000 data points rest all c is still 0
		
		for_n
		//for (c = cent, i = 0; i < n_cluster; i++, c++)
		{ c->x /= c->group; c->y /= c->group; }
 
		
		//problem ends
		
		
		changed = 0;
		/* find closest centroid of each point */
		for_len
		//for (j = 0, p = pts; j < len; j++, p++)
		{
			min_i = nearest(p, cent, n_cluster, 0);
			if (min_i != p->group) {
				changed++;
				p->group = min_i;
			}
		}
	} while (changed > (len >> 5)); /* stop when 99.9% of points are good */
	// now actual groups with neareast centroid is assigned to all data points
	for_n 
	//for (c = cent, i = 0; i < n_cluster; i++, c++)
	{ c->group = i; }
 
	return cent; // from diagram
}
 /*
void print_eps(point pts, int len, point cent, int n_cluster)
{
#	define W 400
#	define H 400
	int i, j;
	point p, c;
	double min_x, max_x, min_y, max_y, scale, cx, cy;
	double *colors = malloc(sizeof(double) * n_cluster * 3);
 
	for_n {
		colors[3*i + 0] = (3 * (i + 1) % 11)/11.;
		colors[3*i + 1] = (7 * i % 11)/11.;
		colors[3*i + 2] = (9 * i % 11)/11.;
	}
 
	max_x = max_y = -(min_x = min_y = HUGE_VAL);
	for_len {
		if (max_x < p->x) max_x = p->x;
		if (min_x > p->x) min_x = p->x;
		if (max_y < p->y) max_y = p->y;
		if (min_y > p->y) min_y = p->y;
	}
	scale = W / (max_x - min_x);
	if (scale > H / (max_y - min_y)) scale = H / (max_y - min_y);
	cx = (max_x + min_x) / 2;
	cy = (max_y + min_y) / 2;
 
	printf("%%!PS-Adobe-3.0\n%%%%BoundingBox: -5 -5 %d %d\n", W + 10, H + 10);
	printf( "/l {rlineto} def /m {rmoveto} def\n"
		"/c { .25 sub exch .25 sub exch .5 0 360 arc fill } def\n"
		"/s { moveto -2 0 m 2 2 l 2 -2 l -2 -2 l closepath "
		"	gsave 1 setgray fill grestore gsave 3 setlinewidth"
		" 1 setgray stroke grestore 0 setgray stroke }def\n"
	);
	for_n {
		printf("%g %g %g setrgbcolor\n",
			colors[3*i], colors[3*i + 1], colors[3*i + 2]);
		for_len {
			if (p->group != i) continue;
			printf("%.3f %.3f c\n",
				(p->x - cx) * scale + W / 2,
				(p->y - cy) * scale + H / 2);
		}
		printf("\n0 setgray %g %g s\n",
			(c->x - cx) * scale + W / 2,
			(c->y - cy) * scale + H / 2);
	}
	printf("\n%%%%EOF");
	free(colors);
#	undef for_n
#	undef for_len
}*/
#define PTS 10000000// no. of data-points
#define K 12//no. of centroids
#define N 2

void createCsv(point p, char *filename)
{
	//Unhandled exception at 0x00f13007 in K_means_trial.exe: 0xC00000FD: Stack overflow.
int i,j,u;
FILE *fp;
double *c_ptr = (double *)malloc(sizeof(double) * 40000);
double c_arr[40000];

double *d_ptr = (double *)malloc(sizeof(double) * 40000);
double d_arr[40000];

double *e_ptr = (double *)malloc(sizeof(double) * 40000);
double e_arr[40000];

c_ptr=c_arr;
d_ptr=d_arr;
e_ptr=e_arr;

filename = strcat(filename, ".csv");
printf("\n Creating %s.csv file",filename);
fp=fopen(filename, "w+");
//fprintf(fp,"X co-ordinates of all the points");
//fprintf(fp,"Y co-ordinates of all the points");
//fprintf(fp, "\n");
for(i=0; i< 40000;i++)
{
	c_arr[i]= p[i].x; 
	d_arr[i]= p[i].y; 
	e_arr[i]= p[i].group;
}
for (j = 0; j < 40000; j++)
{	
    fprintf(fp,",\n%f ",c_ptr[j]);
	fprintf(fp,",%f ",d_ptr[j]);
	fprintf(fp,",%f ",e_ptr[j]);
}
fclose(fp);
printf("\n %sfile created",filename);
}
void grouped_points(point p, char *filename)
{
	int i,j,u;
	FILE *fp;
	double *c_ptr = (double *)malloc(sizeof(double) * K);
	double c_arr[12];
	c_ptr=c_arr;
	//char *filename=str;
	printf("\n Creating %s.csv file",filename);
	
	//str=strcat(str,".csv");
	filename = strcat(filename, ".csv");
	fp=fopen(filename,"w+");
	fprintf(fp,"No. of Points in Cluster");
	for(i=0; i < K; i++)
	{	c_arr[i]=0; 
	}
		for (j = 0; j < PTS; j++, p++)
		{
			u=p->group;
			c_arr[u]+=1;
		}
		for (j = 0; j < K; j++)
	{	
    fprintf(fp,",\n%f ",c_ptr[j]);
	}
fclose(fp);
printf("\n %sfile created",filename);
}
//point v = gen_xy(PTS, 10);
int main()
{
	clock_t start, end;
	double cpu_time_used;
	int i, j;
	char str[100]= "CSVForBarPlot";
	char str1[100] = "CSVForScatterPlot";
	point v,c;
	start = clock();
	v = gen_xy(PTS, 10);//v is the pointer which points to the first element of the data-set
	c = lloyd(v, PTS, K);

	grouped_points(v, str);
	createCsv(v,str1);
	end= clock();
	printf("\nC code successfully completed.");
	cpu_time_used= ((double)(end-start))/CLOCKS_PER_SEC;
	printf("\n\n\nCpu time used is : %f", cpu_time_used);
	getch();
	
	return 0;
	
} 
