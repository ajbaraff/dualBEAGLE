#include "queue.h"

struct queue *new_queue(void)
{
  struct queue *new = malloc(sizeof(struct queue));
  new->front = new->rear = NULL;

  return new;
}

void free_queue(struct queue *queue)
{
  while(!empty(queue)) dequeue(queue);
  free(queue);
}

void enqueue(struct queue *queue, void *node)
{
  struct queue_node *new = malloc(sizeof(struct queue_node));
  new->node = node;
  new->next = NULL;
  if(empty(queue)) queue->front = queue->rear = new;
  else {
    queue->rear->next = new;
    queue->rear = new;
  }
}

void *dequeue(struct queue *queue)
{
  struct queue_node *temp = queue->front;
  void *ret;
  if(empty(queue)) return NULL;
  else {
    queue->front = queue->front->next;
    ret = temp->node;
    free(temp);
    return ret;
  }
}

int empty(struct queue *queue)
{
  return queue->front == NULL;
}
