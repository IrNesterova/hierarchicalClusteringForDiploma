#include "postgres.h"

#include "executor/spi.h"
#include "utils/builtins.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "funcapi.h"
#include <stdint.h>

#define NOT_USED  0 /* node is currently not used */
#define LEAF_NODE 1 /* node contains a leaf node */
#define A_MERGER  2 /* node contains a merged pair of root clusters */
#define A_DIVISIVE 3 /* node contains a split pair of clusters */
#define MAX_LABEL_LEN 16

#define AVERAGE_LINKAGE  1 /* choose average distance */
#define CENTROID_LINKAGE 2 /* choose distance between cluster centroids */
#define COMPLETE_LINKAGE 3 /* choose maximum distance */
#define SINGLE_LINKAGE   4 /* choose minimum distance */

#define alloc_mem(N, T) (T *) calloc(N, sizeof(T))
#define alloc_fail(M) elog(ERROR, "Failed to allocate memory for %s. \n", M)
#define invalid_node(I) elog(ERROR, "Invalid cluster node at index %d. \n", I)

#define EUCLEDIAN_DISTANCE 1
#define CHEBYSHEV_DISTANCE 2
#define MANHATTAN_DISTANCE 3
#define COSINE_SIMILARITY 4

#define AGNES 'a'
#define DIANA 'd'

typedef struct cluster_s cluster_t;
typedef struct cluster_node_s cluster_node_t;
typedef struct neighbour_s neighbour_t;
typedef struct item_s item_t;

typedef struct cluster_node_a_s cluster_node_a_t;
typedef struct cluster_node_d_s cluster_node_d_t;
int clusterN = 1;

char *tableName;
char *commandsArray = "select * from hier_test_2; select * from hier_test;";

float (*distance_fptr)(float **, const int *, const int *, int, int);


typedef struct array_for_hclust{
	char label[MAX_LABEL_LEN];
	int id;
} harray;
typedef struct coord_s {
	float x,y;
	float coord_array[100];
        char *char_array[100];
} coord_t;

float (*distance_method) (const coord_t*, const coord_t*);

typedef struct size_t size;
typedef struct Real real;


struct cluster_s {
        int num_items; /* number of items that was clustered */
        int num_clusters; /* current number of root clusters */
        int num_nodes; /* number of leaf and merged clusters */
        cluster_node_a_t *nodes; /* leaf and merged clusters */
	cluster_node_d_t *nodesD; 
        float **distances; /* distance between leaves */
};


struct cluster_node_s {
	cluster_node_a_t *nodes_a;
	cluster_node_d_t *nodes_d;
        int type; /* type of the cluster node */
        int is_root; /* true if cluster hasn't merged with another */
        int height; /* height of node from the bottom */
        coord_t centroid; /* centroid of this cluster */
        char *label; /* label of a leaf node */
        int *merged; /* indexes of root clusters merged */
        int num_items; /* number of leaf nodes inside new cluster */
        int *items; /* array of leaf nodes indices inside merged clusters */
        neighbour_t *neighbours; /* sorted linked list of distances to roots */
};

struct cluster_node_a_s{
	int type;
	int is_root;
	int height;
	coord_t centroid;
	char *label;
	int *merged;
	int num_items;
	int *items;
	neighbour_t *neighbours;
};

struct cluster_node_d_s{
 	item_t* items;
	int num_items;
	char *label;
	coord_t centroid;
	int is_root;
	int type;
	int *merged;
};

struct neighbour_s {
        int target; /* the index of cluster node representing neighbour */
        float distance; /* distance between the nodes */
        neighbour_t *next, *prev; /* linked list entries */
};

struct item_s {
        coord_t coord; /* coordinate of the input data point */
        char label[MAX_LABEL_LEN]; /* label of the input data point */
};

cluster_t* (*hierarchy) (int, item_t*);
float euclidean_distance(const coord_t *a, const coord_t *b)
{	
        float withoutsqrt = 0;
        for (int i = 0; i < 100; i++){
                withoutsqrt += pow(a->coord_array[i] - b->coord_array[i], 2);
        }
        float smth = sqrt(withoutsqrt);
        return smth;	
}

float jaccard(const coord_t *a, const coord_t *b)

{
        int in = 0;
        int un = 0;
        for (int i = 0; i < 100; i++){
            uint8_t bitmap1[256] = {0};
            uint8_t bitmap2[256] = {0};
            for (int j = 0; j < 256; j++){
                    bitmap1[a->char_array[i][j]] = 1;
            }    
            for (int j = 0; j < 256; j++){
                    bitmap2[b->char_array[i][j]] = 1;
            }
            for (int i = 0; i < 256; i++){
                    in += bitmap1[i]&&bitmap2[i];
                    un+= bitmap1[i] || bitmap2[i];
            }
        }
        float jacc = (float)in/un*100.0;
        return jacc;
}

float soresen_coeff(const coord_t *a, const coord_t *b){
        float coeffsum = 0;
        for (int i = 0; i < 100; i++){
                size_t length1 = strlen(a->char_array[i]) - 1;
                size_t length2 = strlen(b->char_array[i]) - 1;
                double matches;
                int k = 0;
                int j = 0;
                while (k < length1 && j <length2){
                        char ac[3] = {a->char_array[k], a->char_array[k+1],'\0'};
                        char bc[3] = {b->char_array[k], b->char_array[k+1], '\0'};
                        int cmp = strcasecmp(ac, bc);
                        if (cmp == 0)
                                matches += 2;
                        k++;
                        j++;
                }
                coeffsum += matches / (length1 + length2);
        }
        return coeffsum;
}
float chebyshev_distance(const coord_t *a, const coord_t *b)
{
        float chebyshev = 0;
        for (int i=0; i < 100; i++){
                if (fabs(a->coord_array[i] - b->coord_array[i]) > chebyshev)
                        chebyshev = fabs(a->coord_array[i] - b->coord_array[i]);
        }
		
	return chebyshev;
}

float manhattan_distance(const coord_t *a, const coord_t *b)
{
        float manhattan = 0;
        for (int i = 0; i < 100; i++){
                manhattan += fabs(a->coord_array[i] - b->coord_array[i]);
        }
	return manhattan;
}

//need to do cosine_similarity
float cosine_distance(const coord_t *a, const coord_t *b)
{       float cos_numerator = 0;
        float cos_denominator_a = 0;
        float cos_denominator_b = 0;
        for (int i = 0; i < 100; i++){
                cos_numerator += a->coord_array[i] * b->coord_array[i];
                cos_denominator_a += pow(a->coord_array[i], 2);
                cos_denominator_b += pow(b->coord_array[2], 2);
        }
	cos_denominator_a = sqrt(cos_denominator_a);
        cos_denominator_b = sqrt(cos_denominator_b);
	float cos_denominator = cos_denominator_a * cos_denominator_b;
	if (cos_denominator == 0){
		return 1;
	} else {
		float cosine = cos_numerator / cos_denominator;
		return cosine;
	}	   		
}

int levenstein(const coord_t *a, const coord_t *b)
{
   
   

}

int min(int a, int b, int c)
{	
	if(a <= b && a <= c)
	{
		return a;
	}
	else if(b <= a && b <= c)
	{
		return b;
	}
	else if(c <= a && c <= b)
	{
		return c;
	}
}
void fill_euclidean_distances(float **matrix, int num_items,
                              const item_t items[])
{

        for (int i = 0; i < num_items; ++i)
                for (int j = 0; j < num_items; ++j) {
                        matrix[i][j] =
                                distance_method(&(items[i].coord),
                                                   &(items[j].coord));
                        matrix[j][i] = matrix[i][j];
                }
}

void fill_dissimilarity_matrix(float **matrix, int num_items, const item_t items[])
{
	for (int i = 0; i < num_items; i++)
		for (int j = 0; j < num_items; j++){
			float similar = distance_method(&(items[i].coord), 
								&(items[j].coord)); 
			matrix[i][j] = similar - M_E;
			matrix[j][i] = matrix[i][j];
		}
}

float **generate_distance_matrix(int num_items, const item_t items[])
{
	
        float **matrix = alloc_mem(num_items, float *);
        if (matrix) {
                for (int i = 0; i < num_items; ++i) {
                        matrix[i] = alloc_mem(num_items, float);
                        if (!matrix[i]) {
                                alloc_fail("distance matrix row");
                                num_items = i;
                                for (i = 0; i < num_items; ++i)
                                        free(matrix[i]);
                                free(matrix);
                                matrix = NULL;
                                break;
                        }
                }
                if (matrix)
			
                        fill_euclidean_distances(matrix, num_items, items);
        } else
                alloc_fail("distance matrix");
        return matrix;
}

float single_linkage(float **distances, const int a[],
                     const int b[], int m, int n)
{
        float min = FLT_MAX, d;
        for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j) {
                        d = distances[a[i]][b[j]];
                        if (d < min)
                                min = d;
                }
	
        return min;
}

float complete_linkage(float **distances, const int a[],
                       const int b[], int m, int n)
{
        float d, max = 0.0 /* assuming distances are positive */;
        for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j) {
                        d = distances[a[i]][b[j]];
                        if (d > max)
                                max = d;
                }
        return max;
}

float average_linkage(float **distances, const int a[],
                      const int b[], int m, int n)
{
        float total = 0.0;
        for (int i = 0; i < m; ++i)
                for (int j = 0; j < n; ++j)
                        total += distances[a[i]][b[j]];
        return total / (m * n);
}

float centroid_linkage(float **distances, const int a[],
                       const int b[], int m, int n)
{
        return 0; /* empty function */
}



float get_distance(cluster_t *cluster, int index, int target)
{
	
        if (index < cluster->num_items && target < cluster->num_items)
                return cluster->distances[index][target];
        else {
                cluster_node_a_t *a = &(cluster->nodes[index]);
                cluster_node_a_t *b = &(cluster->nodes[target]);
                if (distance_fptr == centroid_linkage)
                        return distance_method(&(a->centroid),
                                                  &(b->centroid));
                else return distance_fptr(cluster->distances,
                                          a->items, b->items,
                                          a->num_items, b->num_items);
        }
}

void free_neighbours(neighbour_t *node)
{
        neighbour_t *t;
        while (node) {
                t = node->next;
                free(node);
                node = t;
        }
}

void free_cluster_nodes(cluster_t *cluster)
{
        for (int i = 0; i < cluster->num_nodes; ++i) {
                cluster_node_a_t *node = &(cluster->nodes[i]);
                if (node->label)
                        free(node->label);
                if (node->merged)
                        free(node->merged);
                if (node->items)
                        free(node->items);
                if (node->neighbours)
                        free_neighbours(node->neighbours);
        }
        free(cluster->nodes);
}

void free_cluster(cluster_t * cluster)
{
        if (cluster) {
                if (cluster->nodes)
                        free_cluster_nodes(cluster);
                if (cluster->distances) {
                        for (int i = 0; i < cluster->num_items; ++i)
                                free(cluster->distances[i]);
                        free(cluster->distances);
                }
                free(cluster);
        }
}

void insert_before(neighbour_t *current, neighbour_t *neighbours,
                   cluster_node_a_t *node)
{
        neighbours->next = current;
        if (current->prev) {
                current->prev->next = neighbours;
                neighbours->prev = current->prev;
        } else
                node->neighbours = neighbours;
        current->prev = neighbours;
}

void insert_after(neighbour_t *current, neighbour_t *neighbours)
{
        neighbours->prev = current;
        current->next = neighbours;
}

void insert_sorted(cluster_node_a_t *node, neighbour_t *neighbours)
{
	
        neighbour_t *temp = node->neighbours;
        while (temp->next) {
                if (temp->distance >= neighbours->distance) {
                        insert_before(temp, neighbours, node);
                        return;
                }
                temp = temp->next;
        }
        if (neighbours->distance < temp->distance)
                insert_before(temp, neighbours, node);
        else
                insert_after(temp, neighbours);
}


neighbour_t *add_neighbour(cluster_t *cluster, int index, int target)
{
	//elog(INFO, "ADDING NEIGHBOUR:%d to %d", index, target);
        neighbour_t *neighbour = alloc_mem(1, neighbour_t);
        if (neighbour) {
                neighbour->target = target;
                neighbour->distance = get_distance(cluster, index, target);
                cluster_node_a_t *node = &(cluster->nodes[index]);
                if (node->neighbours)
                        insert_sorted(node, neighbour);
                else
                        node->neighbours = neighbour;
        } else
                alloc_fail("neighbour node");
        return neighbour;
}

cluster_t *update_neighbours(cluster_t *cluster, int index)
{
	
        cluster_node_a_t *node = &(cluster->nodes[index]);
        if (node->type == NOT_USED) {
                invalid_node(index);
                cluster = NULL;
        } else {
                int root_clusters_seen = 1, target = index;
                while (root_clusters_seen < cluster->num_clusters) {
                        cluster_node_a_t *temp = &(cluster->nodes[--target]);
                        if (temp->type == NOT_USED) {
                                invalid_node(index);
                                cluster = NULL;
                                break;
                        }
                        if (temp->is_root) {
                                ++root_clusters_seen;
				
                                add_neighbour(cluster, index, target);
                        }
                }
        }
        return cluster;
}

#define init_leaf(cluster, node, item, len)             \
        do {                                            \
                strncpy(node->label, item->label, len); \
                node->centroid = item->coord;           \
                node->type = LEAF_NODE;                 \
                node->is_root = 1;                      \
                node->height = 0;                       \
                node->num_items = 1;                    \
                node->items[0] = cluster->num_nodes++;  \
        } while (0)                                     \

cluster_node_a_t *add_leaf(cluster_t *cluster, const item_t *item)
{
        cluster_node_a_t *leaf = &(cluster->nodes[cluster->num_nodes]);
        int len = strlen(item->label) + 1;
        leaf->label = alloc_mem(len, char);
        if (leaf->label) {
                leaf->items = alloc_mem(1, int);
                if (leaf->items) {
                        init_leaf(cluster, leaf, item, len);
                        cluster->num_clusters++;
                } else {
                        alloc_fail("node items");
                        free(leaf->label);
                        leaf = NULL;
                }
        } else {
                alloc_fail("node label");
                leaf = NULL;
        }
        return leaf;
}

#undef init_leaf

cluster_t *add_leaves(cluster_t *cluster, item_t *items)
{
	
        for (int i = 0; i < cluster->num_items; ++i) {
		
                if (add_leaf(cluster, &items[i]))
                        update_neighbours(cluster, i);
                else {
                        cluster = NULL;
                        break;
                }
        }
        return cluster;
}

void print_cluster_items(cluster_t *cluster, int index)
{

        cluster_node_a_t *node = &(cluster->nodes[index]);
	elog(INFO, "cluster number:%d", clusterN);
        elog(INFO, "Items: ");
        if (node->num_items > 0) {
                char command[200];
		sprintf(command, "update %s set cluster_id = %d where id = %s",tableName, clusterN, cluster->nodes[node->items[0]].label);
 		
		SPI_exec(command, 1);
		elog(INFO, "label: %s", cluster->nodes[node->items[0]].label);
                for (int i = 1; i < node->num_items; ++i){
			sprintf(command, "update %s set cluster_id = %d where id = %s", tableName, clusterN, cluster->nodes[node->items[i]].label);
 		
		SPI_exec(command, 1);
		elog(INFO, ", %s", cluster->nodes[node->items[i]].label);
}
        }	
	clusterN = clusterN + 1;	
}

void print_cluster_node(cluster_t *cluster, int index)
{
        cluster_node_a_t *node = &(cluster->nodes[index]);
        print_cluster_items(cluster, index);
        
        neighbour_t *t = node->neighbours;
        while (t) {
               
                t = t->next;
        }
}

void merge_items(cluster_t *cluster, cluster_node_a_t *node,
                 cluster_node_a_t **to_merge)
{
	//elog(INFO, "MERGING ITEMS");
        node->type = A_MERGER;
        node->is_root = 1;
        node->height = -1;

        /* copy leaf indexes from merged clusters */
        int k = 0, idx;
        
        //calculate centroids

        coord_t centroid = { .coord_array = 0};
        for (int i = 0; i < 2; ++i) {
                cluster_node_a_t *t = to_merge[i];
                t->is_root = 0; /* no longer root: merged */
                if (node->height == -1 ||
                    node->height < t->height)
                        node->height = t->height;
                for (int j = 0; j < t->num_items; ++j) {
                        idx = t->items[j];
                        node->items[k++] = idx;
                }
                for (int j = 0; j < 100; j++)
                        centroid.coord_array[i] += t->num_items * t->centroid.coord_array[i];
                        //elog(NOTICE, "centroid: i:%d, val:%lf", i, centroid.coord_array[i]);
             
        }
        /* calculate centroid */
        for (int j = 0; j < 100; j++){
                node->centroid.coord_array[j] = centroid.coord_array[j] / k;        
        }
        node->height++;
}

#define merge_to_one(cluster, to_merge, node, node_idx)         \
        do {                                                    \
                node->num_items = to_merge[0]->num_items +      \
                        to_merge[1]->num_items;                 \
                node->items = alloc_mem(node->num_items, int);  \
                if (node->items) {                              \
                        merge_items(cluster, node, to_merge);   \
                        cluster->num_nodes++;                   \
                        cluster->num_clusters--;                \
                        update_neighbours(cluster, node_idx);   \
                } else {                                        \
                        alloc_fail("array of merged items");    \
                        free(node->merged);                     \
                        node = NULL;                            \
                }                                               \
        } while(0)                                              \

cluster_node_a_t *merge(cluster_t *cluster, int first, int second)
{
        int new_idx = cluster->num_nodes;
        cluster_node_a_t *node = &(cluster->nodes[new_idx]);
        node->merged = alloc_mem(2, int);
        if (node->merged) {
                cluster_node_a_t *to_merge[2] = {
                        &(cluster->nodes[first]),
                        &(cluster->nodes[second])
                };
                node->merged[0] = first;
                node->merged[1] = second;
                merge_to_one(cluster, to_merge, node, new_idx);
        } else {
                alloc_fail("array of merged nodes");
                node = NULL;
        }
        return node;
}

#undef merge_to_one

void find_best_distance_neighbour(cluster_node_a_t *nodes,
                                  int node_idx,
                                  neighbour_t *neighbour,
                                  float *best_distance,
                                  int *first, int *second)
{
        while (neighbour) {
                if (nodes[neighbour->target].is_root) {
                        if (*first == -1 ||
                            neighbour->distance < *best_distance) {
                                *first = node_idx;
                                *second = neighbour->target;
                                *best_distance = neighbour->distance;
                        }
                        break;
                }
                neighbour = neighbour->next;
        }
}


int find_clusters_to_merge(cluster_t *cluster, int *first, int *second)
{
        float best_distance = 0.0;
        int root_clusters_seen = 0;
        int j = cluster->num_nodes; /* traverse hierarchy top-down */
        *first = -1;
        while (root_clusters_seen < cluster->num_clusters) {
                cluster_node_a_t *node = &(cluster->nodes[--j]);
                if (node->type == NOT_USED || !node->is_root)
                        continue;
                ++root_clusters_seen;
                find_best_distance_neighbour(cluster->nodes, j,
                                             node->neighbours,
                                             &best_distance,
                                             first, second);
        }
        return *first;
}

cluster_t *merge_clusters(cluster_t *cluster)
{
        int first, second;
        while (cluster->num_clusters > 1) {
                if (find_clusters_to_merge(cluster, &first, &second) != -1)
                        merge(cluster, first, second);
        }
        return cluster;
}

#define init_cluster(cluster, num_items, items)                         \
        do {                                                            \
		elog(INFO, "INITIATING CLUSTER");                       \
                                                                        \
                cluster->distances =                                    \
                        generate_distance_matrix(num_items, items);     \
                if (!cluster->distances)                                \
                        goto cleanup;                                   \
                cluster->num_items = num_items;                         \
                cluster->num_nodes = 0;                                 \
                cluster->num_clusters = 0;                              \
		elog(INFO, "trying to add leaves");                     \
                if (add_leaves(cluster, items))                         \
                        merge_clusters(cluster);                        \
                else                                                    \
                        goto cleanup;                                   \
        } while (0)                                                     \


cluster_t *agglomerate(int num_items, item_t *items)
{
        cluster_t *cluster = alloc_mem(1, cluster_t);
	elog(INFO, "ALLOCATING MEMORY TO CLUSTER");
        if (cluster) {
                cluster->nodes = alloc_mem(2 * num_items - 1, cluster_node_a_t);
                if (cluster->nodes)
                        init_cluster(cluster, num_items, items);
                else {
                        alloc_fail("cluster nodes");
                        goto cleanup;
                }
        } else
                alloc_fail("cluster");
        goto done;

cleanup:
        free_cluster(cluster);
        cluster = NULL;

done:
        return cluster;
}
#undef init_cluster

cluster_node_d_t *add_item(cluster_t *cluster, item_t *items){
        
	size_t lengthItem = sizeof(&(items))/sizeof(&(items[0]));
	
	
        for (int i = 0; i < lengthItem; i++){
		elog(NOTICE, "TRYING TO ADD ITEMS TO DIANA");
                cluster->nodesD[0].items[i] = items[i];
		elog(NOTICE, "ADDED ONE ITEM");
        }
        cluster->nodesD[0].num_items = lengthItem;
        
        elog(NOTICE, "pls help, there are %d items:", cluster->nodesD[0].num_items); 
}
cluster_t *add_into_cluster(cluster_t *cluster, item_t *items)
{
	elog(NOTICE, "ADDING INTO CLUSTERS");
	add_item(cluster, items);	
}

#define init_cluster_decisive(cluster, num_items, items)\
        cluster->num_clusters = 1;                      \
        cluster->num_items = num_items;                 \
	cluster->num_nodes = 1;				\
	add_into_cluster(cluster, items);		\
	do {						\
                               \
	} while (0)					\
	

cluster_t *decisive(int num_items, item_t *items)
{
        elog(NOTICE, "DECISIVE ALGORITHM");
	cluster_t *cluster = alloc_mem(1, cluster_t);
	cluster->nodesD = alloc_mem(2*num_items - 1, cluster_node_d_t);
	if (cluster->nodesD){
		init_cluster_decisive(cluster, num_items, items);
	}else {	
		alloc_fail("cluster nodes");
		goto cleanup;
	}
	
	goto done;

cleanup:
	free_cluster(cluster);
	cluster = NULL;
done:
	return cluster;
}
	//do init cluster with 1 node, and all items in it

#undef init_cluster_decisive

int print_root_children(cluster_t *cluster, int i, int nodes_to_discard)
{
        cluster_node_a_t *node = &(cluster->nodes[i]);
        int roots_found = 0;
        if (node->type == A_MERGER) {
                for (int j = 0; j < 2; ++j) {
			
                        int t = node->merged[j];
                        if (t < nodes_to_discard) {
                                print_cluster_items(cluster, t);
                                ++roots_found;
                        }
                }
        }
        return roots_found;
}

void get_k_clusters(cluster_t *cluster, int k)
{
        int i = cluster->num_nodes - 1;
        int roots_found = 0;
        int nodes_to_discard = cluster->num_nodes - k + 1;
	elog(INFO, "cluster->num-nodes:i:%d", i);
        while (k) {
		
                if (i < nodes_to_discard) {	
                        print_cluster_items(cluster, i);
                        roots_found = 1;
		
                } else
                        roots_found = print_root_children(cluster, i,
                                                          nodes_to_discard);
                k -= roots_found;
      
                --i;
        }
}

void print_cluster(cluster_t *cluster)
{
        for (int i = 0; i < cluster->num_nodes; ++i)
                print_cluster_node(cluster, i);
}

int read_items(item_t *items, float coord_arr[], int call, int arrsize)
{	
       
	item_t *t = &(items[call]); 
        char label1[16];


	int label = (int) coord_arr[0];
	elog(NOTICE, "%d", label);

        sprintf(label1, "%d", label);
      
	sscanf(label1, "%s", t->label);
	for (int i = 1; i < arrsize; i++){
		t->coord.coord_array[i] = coord_arr[i];                
	}
	
	
        return -1;
}

void set_linkage(int linkage_type)
{
        switch (linkage_type) {
        case AVERAGE_LINKAGE:
                distance_fptr = average_linkage;
		elog(INFO, "setting linkage as average linkage");
                break;
        case COMPLETE_LINKAGE:
                distance_fptr = complete_linkage;
		elog(INFO, "setting linkage as complete linkage");
                break;
        case CENTROID_LINKAGE:
                distance_fptr = centroid_linkage;
		elog(INFO, "setting linkage as centroid linkage");
                break;
        case SINGLE_LINKAGE:
        default: distance_fptr = single_linkage;
	elog(INFO, "setting as single linkage");
        }
}

void set_distance(int distance)
{
	switch (distance){
		case EUCLEDIAN_DISTANCE:
			distance_method = euclidean_distance;
			elog(INFO, "setting distance method as eucledian");
			break;
		case CHEBYSHEV_DISTANCE:
			distance_method = chebyshev_distance;
			elog(INFO, "setting distance method as chebyshev");
			break;
		case MANHATTAN_DISTANCE:
			distance_method = manhattan_distance;
			elog(INFO, "setting distance method as manhattan");
			break;
		case COSINE_SIMILARITY:
			distance_method = cosine_distance;
			elog(INFO, "setting distance as cosine similarity");
			break;
		default:
			distance_method = euclidean_distance;
			break;
	}

}			

void set_hierarchy(char* method){
	//elog(ERROR, "hierarchy method:%s, %d", method, strlen(method));
	switch (method[0]) {
		case 97 :
			hierarchy = agglomerate;
			//matrix = fill_euclidean_distances;
			break;
		case 65:
			hierarchy = agglomerate;
			//matrix = fill_euclidean_distances;
			break;
		case 100:
			hierarchy = decisive;
			//matrix = fill_dissimilarity_matrix;
			break;
		case 68:
			hierarchy = decisive;
			//matrix = fill_dissimilarity_matrix;
			break;
		default:
			hierarchy = agglomerate;
			//matrix = fill_euclidean_distances;
			break;
	}

}
	
int process_input(item_t **items, float row_arr[], int proc, int j, int arrsize)
	
{
	if (j == 0) {
        	*items = alloc_mem(proc, item_t);
		
	}
	read_items(*items, row_arr, j, arrsize);
        return proc;
	
}

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(hclust);

Datum
hclust(PG_FUNCTION_ARGS)
{
    //methods for linkage and distance
    int linkageMethod;
    int distanceMethod;

    char *command;
    int cnt;
    int ret;
    int ret1;
    uint64 proc;

    //int for cluster number
    int clusternumber;
	tableName = text_to_cstring(PG_GETARG_TEXT_PP(2));
    /* Convert given text object to a C string */
    command = text_to_cstring(PG_GETARG_TEXT_PP(0));
    
    SPI_connect();
    char getMaxAndMin[200];

    char commandforaltertable[200];
    sprintf(commandforaltertable, "alter table %s add column if not exists cluster_id integer", tableName);
    SPI_exec(commandforaltertable, 1);
    ret = SPI_exec(command, 0);

    proc = SPI_processed;
    //cluster number
    clusternumber = PG_GETARG_INT32(1);
    
    int spi;
    item_t *items = NULL;	
    int num_items = proc;
  
    linkageMethod = PG_GETARG_INT32(3);
    char table_label[MAX_LABEL_LEN];   
 
    distanceMethod = PG_GETARG_INT32(4);
    char* clusteringType = text_to_cstring(PG_GETARG_TEXT_PP(5));
    set_hierarchy(clusteringType);
    set_linkage(linkageMethod);
    set_distance(distanceMethod);
    int call_cntr;
    int max_calls;

    float row_arr[100] = {'\0'};

    clusterN = 1;
    if (clusternumber < 1) 
	elog(ERROR, "Cluster number < 1");
    if (ret > 0 && SPI_tuptable != NULL)
    {
        TupleDesc tupdesc = SPI_tuptable->tupdesc;
        SPITupleTable *tuptable = SPI_tuptable;
       
        uint64 j;
	char* type;
	
	char *columnName;
        for (j = 0; j < proc; j++)
        {
            HeapTuple tuple = tuptable->vals[j];
            int i;    
	    for (i = 1; i <= tupdesc->natts; i++){
		

		 char* columnType = SPI_getvalue(tuple, tupdesc, i);	
			
		 char* value;
		 value = SPI_getvalue(tuple, tupdesc, i);
		
		 if (value == NULL){
			elog(NOTICE, "value of %d row of %i column is null, replacing with 0");
			row_arr[i-1] = 0;
                       
		} else {
                	row_arr[i-1] = atof(value);
		}
		
		                
	    }
		process_input(&items, row_arr, proc, j, tupdesc->natts); 
	}	  	
	elog(NOTICE, "Input processed");
        
        cluster_t *cluster = hierarchy(num_items, items);
	elog(INFO, "num items: %d", cluster->num_items);

	if (cluster)
			
		elog(INFO, "\n\n%d CLUSTERS\n"
		                                "--------------------\n", clusternumber);
		get_k_clusters(cluster, clusternumber);			
        }
	
	SPI_finish();	
    PG_RETURN_INT64(proc);

}



