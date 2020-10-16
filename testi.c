main()
{
int i, *p;
char c;
short s;
long l;
p=&i;
printf("character size is %d\n",sizeof(c));
printf("short     size is %d\n",sizeof(s));
printf("integer   size is %d\n",sizeof(i));
printf("long      size is %d\n",sizeof(l));
printf("pointer   size is %d\n",sizeof(p));
printf("pointer  value is %lx\n",(long)p);
}
