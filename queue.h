#ifndef queue_H_
#define queue_H_

#include <stdlib.h>

struct queue_node
{
  void *node;
  struct queue_node *next;
};

struct queue
{
  struct queue_node *front;
  struct queue_node *rear;
};

struct queue *new_queue(void);
void enqueue(struct queue*, void*);
void *dequeue(struct queue*);
int empty(struct queue*);

#endif // queue_H_
