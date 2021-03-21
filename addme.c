#include "postgres.h"

#include "executor/spi.h"
#include "utils/builtins.h"
#include <string.h>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <stdio.h>
#include "funcapi.h"

#define NOT_USED  0 /* node is currently not used */
#define LEAF_NODE 1 /* node contains a leaf node */
#define A_MERGER  2 /* node contains a merged pair of root clusters */
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

int clusterN = 1;

float (*distance_fptr)(float **, const int *, const int *, int, int);



typedef struct coord_s {
        float x, y, z, w, v, u, t, r, s, q;
} coord_t;

float (*distance_method) (const coord_t*, const coord_t*);



struct cluster_s {
        int num_items; /* number of items that was clustered */
        int num_clusters; /* current number of root clusters */
        int num_nodes; /* number of leaf and merged clusters */
        cluster_node_t *nodes; /* leaf and merged clusters */
        float **distances; /* distance between leaves */
};

struct cluster_node_s {
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
        float smth = sqrt(pow(a->x - b->x, 2) + pow(a->y - b->y, 2) + pow(a->z - b->z, 2) + pow(a->w - b->w, 2) + pow(a->v - b->v, 2) + pow(a->u - b->u, 2) + pow(a->t - b->t, 2) + pow (a->r - b->r, 2) + pow(a->s - b->s, 2) + pow(a->q - b->q,2));
	return smth;	
}

float chebyshev_distance(const coord_t *a, const coord_t *b)
{
	float chebyshev = fmaxf(fabs(a->x - b->x), fmaxf(fabs(a->y - b->y), fmaxf(fabs(a->z - b->z), fmaxf(fabs(a->w - b->w), fmaxf(fabs(a->v - b->v), fmaxf(fabs(a->u - b->u), fmaxf(fabs(a->t - b->t), fmax(fabs(a->r - b->r), fmax(fabs(a->s - b->s), fabs(a->q - b->q))))))))));
	elog(INFO, "chebyshev distance:%lf", chebyshev);	
	return chebyshev;
}

float manhattan_distance(const coord_t *a, const coord_t *b)
{
	float manhattan =  fabs(a->x - b->x) + fabs(a->y - b->y) + fabs(a->z - b->z) + fabs(a->w - b->w) + fabs(a->u - b->u) + fabs(a->v - b->u) + fabs(a->t - b->t) + fabs(a->r - b->r) + fabs(a->s - b->s) + fabs(a->q - b->q);
	elog(INFO, "manhattan:%lf", manhattan);
	return manhattan;
}

//need to do cosine_similarity
float cosine_distance(const coord_t *a, const coord_t *b)
{
	float cos_numerator = (a->x * b->x + a->y * b->y + a->z * b->z + a->w * b->w + a->u * b->u + a->v * b->v + a->t * b->t + a->r * b->r + a->s * b->s + a->q * b->q);
	float cos_denominator_a = sqrt(pow(a->x,2) + pow(a->y, 2) + pow(a->z, 2) + pow(a->w, 2) + pow(a->u,2) + pow(a->v, 2) + pow(a->t, 2) + pow (a->s, 2) + pow(a->q,2) + pow(a->r,2));
	float cos_denominator_b = sqrt(pow(b->x,2) + pow(b->y, 2) + pow(b->z, 2) + pow(b->w, 2) + pow(b->u,2) + pow(b->v, 2) + pow(b->t, 2) + pow (b->s, 2) + pow(b->q,2) + pow(b->r,2)); 
	float cos_denominator = cos_denominator_a * cos_denominator_b;
	if (cos_denominator == 0){
		return 1;
	} else {
		float cosine = cos_numerator / cos_denominator;
		return cosine;
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
	for (int i = 0; i < num_items; i++){
		elog(INFO, "%s", items[i].label);
	} 
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
	elog(INFO, "%f", min);
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
                cluster_node_t *a = &(cluster->nodes[index]);
                cluster_node_t *b = &(cluster->nodes[target]);
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
                cluster_node_t *node = &(cluster->nodes[i]);
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
                   cluster_node_t *node)
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

void insert_sorted(cluster_node_t *node, neighbour_t *neighbours)
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
	elog(INFO, "ADDING NEIGHBOUR:%d to %d", index, target);
        neighbour_t *neighbour = alloc_mem(1, neighbour_t);
        if (neighbour) {
                neighbour->target = target;
                neighbour->distance = get_distance(cluster, index, target);
                cluster_node_t *node = &(cluster->nodes[index]);
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
	
        cluster_node_t *node = &(cluster->nodes[index]);
        if (node->type == NOT_USED) {
                invalid_node(index);
                cluster = NULL;
        } else {
                int root_clusters_seen = 1, target = index;
                while (root_clusters_seen < cluster->num_clusters) {
                        cluster_node_t *temp = &(cluster->nodes[--target]);
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

cluster_node_t *add_leaf(cluster_t *cluster, const item_t *item)
{
        cluster_node_t *leaf = &(cluster->nodes[cluster->num_nodes]);
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
        cluster_node_t *node = &(cluster->nodes[index]);
	elog(INFO, "global parameter:%d", clusterN);
        elog(INFO, "Items: ");

        if (node->num_items > 0) {
		
		elog(INFO, "%s", cluster->nodes[node->items[0]].label);
		
		
                for (int i = 1; i < node->num_items; ++i)
			elog(INFO, ", %s", cluster->nodes[node->items[i]].label);
                        
        }	
	clusterN = clusterN + 1;
	
}

void print_cluster_node(cluster_t *cluster, int index)
{
        cluster_node_t *node = &(cluster->nodes[index]);
	elog(INFO, "Node %d - height: %d, centroid: (%5.3f, %5.3f)\n",
                index, node->height, node->centroid.x, node->centroid.y);
        if (node->label)
		elog(INFO, "\tLeaf: %s\n\t", node->label);
                
        else
		 elog(INFO, "\tMerged: %d, %d\n\t",
                        node->merged[0], node->merged[1]);
             
        print_cluster_items(cluster, index);
        elog(INFO, "\tNeighbours: ");
        neighbour_t *t = node->neighbours;
        while (t) {
                elog(INFO, "target \n\t\t%2d: and distance %5.3f", t->target, t->distance);
                t = t->next;
        }
}

void merge_items(cluster_t *cluster, cluster_node_t *node,
                 cluster_node_t **to_merge)
{
	//elog(INFO, "MERGING ITEMS");
        node->type = A_MERGER;
        node->is_root = 1;
        node->height = -1;

        /* copy leaf indexes from merged clusters */
        int k = 0, idx;
        coord_t centroid = { .x = 0.0, .y = 0.0 };
        for (int i = 0; i < 2; ++i) {
                cluster_node_t *t = to_merge[i];
                t->is_root = 0; /* no longer root: merged */
                if (node->height == -1 ||
                    node->height < t->height)
                        node->height = t->height;
                for (int j = 0; j < t->num_items; ++j) {
                        idx = t->items[j];
                        node->items[k++] = idx;
                }
                centroid.x += t->num_items * t->centroid.x;
                centroid.y += t->num_items * t->centroid.y;
        }
        /* calculate centroid */
        node->centroid.x = centroid.x / k;
        node->centroid.y = centroid.y / k;
        node->height++;
}

#define merge_to_one(cluster, to_merge, node, node_idx)         \
        do {                                                    \
		elog(INFO, "MERGING CLUSTERS TO ONE"); \
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

cluster_node_t *merge(cluster_t *cluster, int first, int second)
{
        int new_idx = cluster->num_nodes;
        cluster_node_t *node = &(cluster->nodes[new_idx]);
        node->merged = alloc_mem(2, int);
        if (node->merged) {
                cluster_node_t *to_merge[2] = {
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

void find_best_distance_neighbour(cluster_node_t *nodes,
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
                cluster_node_t *node = &(cluster->nodes[--j]);
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
				                                        \
			                                                \
                if (!cluster->distances)                                \
                        goto cleanup;                                   \
                cluster->num_items = num_items;                         \
                cluster->num_nodes = 0;                                 \
                cluster->num_clusters = 0;                              \
		elog(INFO, "trying to add leaves"); \
                if (add_leaves(cluster, items))                 \
                        merge_clusters(cluster);                        \
                else                                                    \
                        goto cleanup;                                   \
        } while (0)                                                     \


cluster_t *agglomerate(int num_items, item_t *items)
{
	
        cluster_t *cluster = alloc_mem(1, cluster_t);
	elog(INFO, "ALLOCATING MEMORY TO CLUSTER");
        if (cluster) {
                cluster->nodes = alloc_mem(2 * num_items - 1, cluster_node_t);
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

cluster_t *add_item(cluster_t *cluster, item_t* item)
{
	
}
cluster_t *add_into_cluster(cluster_t *cluster, item_t *items)
{
	for (int i = 0; i < cluster->num_items; i++){
		add_item(cluster, &(items[i]));
	}
}


#define init_cluster_decisive(cluster, num_items, items) \
	add_into_cluster(cluster, items);		\
	

cluster_t *decisive(int num_items, item_t *items)
{
	cluster_t *cluster = alloc_mem(1, cluster_t);
	cluster->nodes = alloc_mem(2*num_items - 1, cluster_node_t);
	if (cluster->nodes){
		//init_cluster_decisive(cluster, num_items, items);
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
        cluster_node_t *node = &(cluster->nodes[i]);
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
	elog(INFO, "i:%d", i);
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

int read_items(item_t *items, char arr[], int call, int arrsize)
{	
	elog(INFO, "trying to read items");
	item_t *t = &(items[call]);
	switch (arrsize){
		case 1:
			elog(ERROR, "1st column is label");
			break;
		case 2:
			sscanf(arr, "%s%f", t->label, &(t->coord.x));
			
			break;
		case 3:
			sscanf(arr, "%s%f%f", t->label, &(t->coord.x), &(t->coord.y));
			break;
		case 4:
			sscanf(arr, "%s%f%f%f", t->label, &(t->coord.x), &(t->coord.y), &(t->coord.z));
			
			break;
		case 5:
			sscanf(arr, "%s%f%f%f%f", t->label, &(t->coord.x), &(t->coord.y), &(t->coord.z), &(t->coord.w));
			break;
		case 6:
			sscanf(arr, "%s%f%f%f%f%f", t->label, &(t->coord.x), &(t->coord.y), &(t->coord.z), &(t->coord.w), &(t->coord.v));
			break;
		case 7:
			sscanf(arr, "%s%f%f%f%f%f%f", t->label, &(t->coord.x), &(t->coord.y), &(t->coord.z), &(t->coord.w), &(t->coord.v), &(t->coord.u));
			break;
		case 8:
			sscanf(arr, "%s%f%f%f%f%f%f%f", t->label, &(t->coord.x), &(t->coord.y), &(t->coord.z), &(t->coord.w), &(t->coord.v), &(t->coord.u), &(t->coord.t));
			break;
		case 9:
			sscanf(arr, "%s%f%f%f%f%f%f%f%f", t->label, &(t->coord.x), &(t->coord.y), &(t->coord.z), &(t->coord.w), &(t->coord.v), &(t->coord.u), &(t->coord.t), &(t->coord.r));
			break; 
		case 10:
			sscanf(arr, "%s%f%f%f%f%f%f%f%f%f", t->label, &(t->coord.x), &(t->coord.y), &(t->coord.z), &(t->coord.w), &(t->coord.v), &(t->coord.u), &(t->coord.t), &(t->coord.r), &(t->coord.s));
			break;
		case 11:
			sscanf(arr, "%s%f%f%f%f%f%f%f%f%f%f", t->label, &(t->coord.x), &(t->coord.y), &(t->coord.z), &(t->coord.w), &(t->coord.v), &(t->coord.u), &(t->coord.t), &(t->coord.r), &(t->coord.s), &(t->coord.q));
			break;
		default:
			elog (ERROR, "too many parameters");
	}
	elog(INFO, "READING ITEMS IS FINISHED");
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
	

int process_input(item_t **items, char arr[], int proc, int j, int arrsize)
	
{
	if (j == 0) {
        	*items = alloc_mem(proc, item_t);
	}
	read_items(*items, arr, j, arrsize);
        return proc;
}

PG_MODULE_MAGIC;

PG_FUNCTION_INFO_V1(addme);

Datum
addme(PG_FUNCTION_ARGS)
{
    //methods for linkage and distance
    int linkageMethod;
    int distanceMethod;

    char *command;
    int cnt;
    int ret;
    uint64 proc;

    //int for cluster number
    int clusternumber;
    char* tablename;

    /* Convert given text object to a C string */
    command = text_to_cstring(PG_GETARG_TEXT_PP(0));
    
    SPI_connect();

    ret = SPI_exec(command, 0);

    proc = SPI_processed;
    //cluster number
    clusternumber = PG_GETARG_INT32(1);
    tablename = text_to_cstring(PG_GETARG_TEXT_PP(2));
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

    clusterN = 1;
    if (clusternumber < 1) 
	elog(ERROR, "Cluster number < 1");
    if (ret > 0 && SPI_tuptable != NULL)
    {
        TupleDesc tupdesc = SPI_tuptable->tupdesc;
        SPITupleTable *tuptable = SPI_tuptable;
       
        uint64 j;
	char* type;
	char table_arr[8192] = {'\0'};
	char *columnName;
        for (j = 0; j < proc; j++)
        {
            HeapTuple tuple = tuptable->vals[j];
            int i;
	    if (j == 0){
		for (i = 2; i <= tupdesc->natts; i++){
			columnName = SPI_fname(tupdesc, i);
			elog(INFO, "%s", columnName);
		}
	    }
	    for (i = 1; i <= tupdesc->natts; i++){
		
		char* value;
		value = SPI_getvalue(tuple, tupdesc, i);
		if (value == NULL){
			if (i == tupdesc->natts){
				strcat(table_arr, "0");
			} else {
				strcat(table_arr, "0 ");
			}
		} else {
			if (i == tupdesc->natts){
				strcat(table_arr, value);
			} else {
				strcat(table_arr, strcat(value, " "));
			}
		}
	    }
	process_input(&items, table_arr, proc, j, tupdesc->natts); 
	memset(table_arr, 0, sizeof(table_arr));     	   	
	
      if (j == proc - 1){
	cluster_t *cluster = hierarchy(num_items, items);
	elog(INFO, "num items: %d", cluster->num_items);
	if (cluster)
			
		elog(INFO, "\n\n%d CLUSTERS\n"
		                                "--------------------\n", clusternumber);
		get_k_clusters(cluster, clusternumber);			
	}	
    }
	
}

    free(items);     
    SPI_finish();
    PG_RETURN_INT64(proc);

}



