//============================================================
//  File:       im_queue.cpp
//  Author:     John Gauch
//  Date:       Spring 2009
//============================================================

#include "im_queue.h"

//============================================================
im_queue::im_queue()
{
   // Init queue
   Head = NULL;
   Tail = NULL;
   Length = 0;
}

//============================================================ 
im_queue::~im_queue()
{
   // TBA
}

//============================================================
bool im_queue::IsFull()
{
   return (Length >= MAX_LENGTH);
}

//============================================================
bool im_queue::IsEmpty()
{
   return (Length == 0);
}

//============================================================
void im_queue::WaitNotFull()
{
   while (IsFull())
      sleep(1);
}

//============================================================
void im_queue::WaitNotEmpty()
{
   while (IsEmpty())
      sleep(1);
}

//============================================================
void im_queue::Insert(int number, im_color & im)
{
   // Insert node into linked list
   im_node *New = new im_node;
   if (New != NULL)
   {
      // Copy data
      New->Number = number;
      New->Image.Copy(im);
      New->Prev = Tail;
      New->Next = NULL;

      // Update head and tail
      if (Tail != NULL)
	 Tail->Next = New;
      else
	 Head = New;
      Tail = New;
      Length++;
   }
}

//============================================================
void im_queue::Remove(int &number, im_color & im)
{
   // Remove node into linked list
   if (Head != NULL)
   {
      // Copy data
      number = Head->Number;
      im.Copy(Head->Image);

      // Update head and tail
      im_node *Temp = Head;
      Head = Head->Next;
      if (Head != NULL)
	 Head->Prev = NULL;
      else
	 Tail = NULL;
      delete Temp;
      Length--;
   }
}

//============================================================
int im_queue::GetLength()
{
   return Length;
}
