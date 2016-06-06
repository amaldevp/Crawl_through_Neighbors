/******************************************************************
 Paper Name: Crawl through Neighbors: A Simple Curve Reconstruction Algorithm
 Authors: Amal Dev Parakkat and Ramanathan Muthuganapathy
 Published in: Computer Graphics Forum (Proceedings of Eurographics Symposium on Geometry Processing 2016)
 
 Implemented using: CGAL 4.6
 Operating System used: Ubuntu 14.04
 Date of final modification: 06-06-2016
 
 Don't hesitate to contact Amal Dev P. (adp.upasana@gmail.com) or Ramanathan M. (emry01@gmail.com) for any queries or issues related to this code.
******************************************************************/
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include<string.h>
#include<fstream>
#if defined(__APPLE__)
#include <GLUT/glut.h>
#else
#include <GL/glut.h>
#endif
#define NENDS 2
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Projection_traits_xy_3<K>  Gt;
typedef CGAL::Delaunay_triangulation_2<Gt> Delaunay;
typedef K::Point_3   Point;
GLdouble width, height;
int wd,loop,ends1[NENDS][2];
Point shape[500000][2],input_points[500000];
int vertex_count=0,shape_index=0,minx=9999,miny=9999,maxx=0,maxy=0;
Delaunay dt;
/**********PRIORITY QUEUE FUNCTIONS**********/
struct node
{
    float priority;
    Delaunay::Edge_iterator afi;
    struct node *link;
};
class prio_queue
{
private:
    node *front;
public:
    prio_queue()
    {
        front = NULL;
    }
    void pq_insert(float priority,Delaunay::Edge_iterator item)/*Inserting an item into priority queue*/
    {
        node *tmp, *q;
        tmp = new node;
        tmp->afi = item;
        tmp->priority = priority;
        if (front == NULL || priority < front->priority)
        {
            tmp->link = front;
            front = tmp;
        }
        else
        {
            q = front;
            while (q->link != NULL && q->link->priority <= priority)
                q=q->link;
            tmp->link = q->link;
            q->link = tmp;
        }
    }
    Delaunay::Edge_iterator pq_del() /* Deleting an item from priority queue*/
    {
        node *tmp;
        if(front == NULL)
            std::cout<<"Queue Underflow\n";
        else
        {
            tmp = front;
            Delaunay::Edge_iterator afi_new;
            afi_new=tmp->afi;
            front = front->link;
            free(tmp);
            return afi_new;
        }
    }
    bool pq_empty() /*Checks whether the priority queue is empty or not*/
    {
        if(front==NULL)
            return 1;
        return 0;
    }
};
prio_queue pq;
/**********OPENGL FUNCTIONS**********/
void init(void)
{
    width  = 1280.0;
    height = 800.0;
    ends1[0][0] = (int)(0.25*width);
    ends1[0][1] = (int)(0.75*height);
    ends1[1][0] = (int)(0.75*width);
    ends1[1][1] = (int)(0.25*height);
    return;
}
void reshape(int w, int h)
{
    width = (GLdouble) w;
    height = (GLdouble) h;
    glViewport(0, 0, (GLsizei) width, (GLsizei) height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glRotatef(-180,180,0,1);
    glOrtho(minx-20.0,maxx+20.0,miny-20.0,maxy+20.0, -1.f, 1.f);
    return;
}
void kbd(unsigned char key, int x, int y)
{
    switch((char)key) {
        case 'q':
        case 27:
            glutDestroyWindow(wd);
            exit(0);
        default:
            break;
    }
    return;
}
void drawFilledCircle(GLfloat x, GLfloat y, GLfloat radius) /*to display vertices*/
{
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    int i,triangleAmount = 20;
    GLfloat twicePi = 2.0f * 3.14159265;
    glBegin(GL_TRIANGLE_FAN);
    glVertex2f(x, y);
    for(i = 0; i <= triangleAmount;i++)
        glVertex2f(x + (radius * cos(i *  twicePi / triangleAmount)),y + (radius * sin(i * twicePi / triangleAmount)));
    glEnd();
}
void pointset(void)
{
    glLineWidth(6.0);
    glEnable( GL_LINE_SMOOTH );
    glEnable( GL_POLYGON_SMOOTH );
    glHint( GL_LINE_SMOOTH_HINT, GL_NICEST );
    glHint( GL_POLYGON_SMOOTH_HINT, GL_NICEST );
    for(int i=0;i<shape_index;i++) /*Displaying edges*/
    {
        glColor3f((164.0/255), (131.0/255), (196.0/255));
        glBegin(GL_LINES);
        glVertex2f(shape[i][0].x(),shape[i][0].y());
        glVertex2f(shape[i][1].x(),shape[i][1].y());
        glEnd();
    }
    glColor3f(0,0,0);
    for(int i=0;i<vertex_count;i++) /*Displaying vertices*/
        drawFilledCircle(input_points[i].x(),input_points[i].y(),2);/*change the last parameter to increase the point size*/
}
void display(void)
{
    glClear(GL_COLOR_BUFFER_BIT);
    pointset();
    glFlush();
    return;
}
/**********RECONSTRUCTION FUNCTIONS**********/
float distance(Point a, Point b)
{
    return (float)(sqrt(abs(((a.x()-b.x())*(a.x()-b.x())))+abs(((a.y()-b.y())*(a.y()-b.y())))));
}
void insert_to_shape(Point a,Point b)/*inserting edge to the shape*/
{
    shape[shape_index][0]=a;
    shape[shape_index][1]=b;
    shape_index++;
}
bool notsmallest(Delaunay::Vertex_handle a,Delaunay::Vertex_handle b)/*checking whether an edge is a potential seed edge or not*/
{
    Delaunay::Vertex_handle vh1;
    Delaunay::Vertex_circulator vc=dt.incident_vertices(a),done(vc);
    double min1=99999;
    if (vc != 0) {
        do {/*Finding the nearest vertex from vertex a */
            if(distance(vc->point(),a->point())<min1)
            {
                min1=distance(vc->point(),a->point());
                vh1=vc;
            }
        }while(++vc != done);
        if(vh1->point()==b->point())/*checking whether the nearest vertex of 'a' is 'b'*/
            return 0;/*potential seed edge*/
    }
    Delaunay::Vertex_circulator vc1=dt.incident_vertices(b),done1(vc1);
    min1=99999;
    if (vc1 != 0) {
        do {/*Finding the nearest vertex from vertex b */
            if(distance(vc1->point(),b->point())<min1)
            {
                min1=distance(vc1->point(),b->point());
                vh1=vc1;
            }
        }while(++vc1 != done1);
        if(vh1->point()==a->point())/*checking whether the nearest vertex of 'b' is 'a'*/
            return 0;/*potential seed edge*/
    }
    return 1;/*'a' is not the nearest vertex of 'b' and 'b' is not the nearest to 'a'*/
}
int main(int argc, char **argv)
{
    init();
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGBA);
    glutInitWindowSize(250, 250);
    wd = glutCreateWindow("Crawl through Neighbors");
    std::ifstream ifs(argv[1]);
    if(ifs.fail())
    {
        printf("File does not exists\n");
        exit(1);
    }
    std::istream_iterator<Point> begin(ifs);
    std::istream_iterator<Point> end;
    dt.insert(begin, end);/*Compute Delaunay triangulation*/
    Delaunay::Vertex_iterator vi=dt.vertices_begin();
    do{
        if(vi->point().x()>maxx)/*Finding the maximum and minimum x,y coordinates*/
            maxx=vi->point().x();//////////////////////////////////////////////////
        if(vi->point().y()>maxy)///////////////////////////////////////////////////
            maxy=vi->point().y();//////////////////////////////////////////////////
        if(vi->point().x()<minx)///////////////////////////////////////////////////
            minx=vi->point().x();//////////////////////////////////////////////////
        if(vi->point().y()<miny)///////////////////////////////////////////////////
            miny=vi->point().y();//////////////////////////////////////////////////
        input_points[vertex_count]=vi->point();/*saving vertices information in an array*/
        vertex_count++;/*number of vertices in the file*/
        if(vertex_count>500000)
        {
            printf("Input size is too high, please increase the size of input_points and shape array\n");
            exit(1);
        }
        vi++;
    }while(vi!=dt.vertices_end());
    float min=9999;
    Delaunay::Vertex_handle ver_hand1,ver_hand2;
    Delaunay::Edge_iterator ei=dt.edges_begin();
    dt.infinite_vertex()->point().flag=1;
    int i;
    do{
        Delaunay::Face& f = *(ei->first);
        i = ei->second;
        if(!notsmallest(f.vertex(f.cw(i)),f.vertex(f.ccw(i))))/*check whether the edge is a seed edge*/
        {
            float dist=distance(f.vertex(f.cw(i))->point(),f.vertex(f.ccw(i))->point());/*calculate the edge length*/
            pq.pq_insert(dist,ei);/*inserting edge into priority queue*/
        }
        ei++;
    }while(ei!=dt.edges_end());
    ei=pq.pq_del();/*delete first edge from priority queue*/
    Delaunay::Face& face_var = *(ei->first);
    i = ei->second;
    insert_to_shape(face_var.vertex(face_var.cw(i))->point(),face_var.vertex(face_var.ccw(i))->point());/*inserting 'ei' to shape*/
    face_var.vertex(face_var.cw(i))->point().flag=1;/*marking vertices as visited*/
    face_var.vertex(face_var.ccw(i))->point().flag=1;//////////////////////////////
    ver_hand1=face_var.vertex(face_var.cw(i));
    ver_hand2=face_var.vertex(face_var.ccw(i));
    Delaunay::Vertex_handle vh1,vh2,vht1,vht2;
    int finished=0,flag=0;
    int finished_one_component=0,cou=0;
    Delaunay::Vertex_handle vi1=ver_hand1,vi2=ver_hand2;
    while(finished==0)/*run until no more seed edge can be found*/
    {
        
        if(flag!=0)
        {
            int found=0;
            min=99999;
            Delaunay::Edge_iterator ei1;
            while(found==0)
            {
            ei1=pq.pq_del();
                Delaunay::Face& f = *(ei1->first);
                int i = ei1->second;
                Delaunay::Vertex_handle vs = f.vertex(f.cw(i));
                Delaunay::Vertex_handle vt = f.vertex(f.ccw(i));
                if(vs->point().flag==0&&vt->point().flag==0)
                {
                    min=distance(vs->point(),vt->point());
                    ver_hand1=vs;
                    ver_hand2=vt;
                    found=1;
                }
                if(pq.pq_empty())
                    break;
            }
            if(found==0||notsmallest(ver_hand1,ver_hand2))
            {
                finished=1;
                break;
            }
                insert_to_shape(ver_hand1->point(),ver_hand2->point());
                ver_hand1->point().flag=1;
                ver_hand2->point().flag=1;
        }
        flag=0;
        finished_one_component=0;
        while(finished_one_component==0)/*run until curve cannot be grown further*/
        {
      Delaunay::Vertex_circulator vc=dt.incident_vertices(ver_hand1),done(vc);/*ver_hand1 and ver_hand2 are the extremity vertices*/
            float min1=9999,min2=9999;
            int entered=0;
            if (vc != 0) {
                do {
                    if(distance(vc->point(),ver_hand1->point())<min1&&vc->point().flag==0)/*find the shortest unvisited vertex from ver_hand1*/
                    {
                        min1=distance(vc->point(),ver_hand1->point());
                        vht1=vh1;
                        vh1=vc;
                        entered=1;/*The curve can be grown further*/
                        
                    }
                }while(++vc != done);
            }
            float min11=9999,min12=9999;
            int fl=0;
            Delaunay::Vertex_circulator vct1=dt.incident_vertices(ver_hand1),donet1(vct1);
            Delaunay::Vertex_handle vh11=vct1,vht1=vct1;
            if (vct1 != 0) {
                do {
                    fl++;
                    if(distance(vct1->point(),ver_hand1->point())<min11)/*find the shortest and second shortest vertex from ver_hand1*/
                    {
                        min12=min11;
                        vh11=vht1;
                        min11=distance(vct1->point(),ver_hand1->point());
                        vht1=vct1;
                    }
                    else
                        if(distance(vct1->point(),ver_hand1->point())<min12)
                        {
                            min12=distance(vct1->point(),ver_hand1->point());
                            vh11=vct1;
                        }
                }while(++vct1 != donet1);
                if(((vh11==ver_hand2&&flag>3)))/*if the curve is closable*/
                {
                    insert_to_shape(ver_hand1->point(),ver_hand2->point());
                    flag=1;
                    finished_one_component=1;
                    break;
                }
            }
            min2=9999;
            Delaunay::Vertex_circulator vc11=dt.incident_vertices(ver_hand2),donen1(vc11);
            if (vc11 != 0) {
                do {
                    if(distance(vc11->point(),ver_hand2->point())<min2&&vc11->point().flag==0)/*find the shortest unvisited vertex from ver_hand2*/
                    {
                        min2=distance(vc11->point(),ver_hand2->point());
                        vht2=vh2;
                        vh2=vc11;
                        entered=1;/*The curve can be grown further*/
                    }
                }while(++vc11 != donen1);
            }
            min11=9999;
            min12=9999;
            Delaunay::Vertex_circulator vct2=dt.incident_vertices(ver_hand2),donet2(vct2);
            Delaunay::Vertex_handle vh21=vct2,vht2=vct2;
            if (vct2 != 0) {
                do {
                    if(distance(vct2->point(),ver_hand2->point())<min11) /*find the shortest and second shortest vertex from ver_hand2*/
                    {
                        min12=min11;
                        vh21=vht2;
                        min11=distance(vct2->point(),ver_hand2->point());
                        vh21=vht2;
                        vht2=vct2;
                    }
                    else
                        if(distance(vct2->point(),ver_hand2->point())<min12)
                        {
                            min12=distance(vct2->point(),ver_hand2->point());
                            vh21=vct2;
                        }
                }while(++vct2 != donet2);
                if(((vh21==ver_hand1&&flag>3)))/*if the curve is closable*/
                {
                    insert_to_shape(ver_hand1->point(),ver_hand2->point());
                    flag=1;
                    finished_one_component=1;
                    break;
                }
            }
            if(min1<min2&&entered==1&&min1!=9999)/*grow from extremity point ver_hand1*/
            {
                vh1->point().flag=1;
                flag++;
                insert_to_shape(ver_hand1->point(),vh1->point());
                ver_hand1=vh1;
            }
            else
                if(entered==1&&min2!=9999)/*grow from extremity point ver_hand2*/
                {
                    vh2->point().flag=1;
                    flag++;
                    insert_to_shape(ver_hand2->point(),vh2->point());
                    ver_hand2=vh2;
                }
                else
                {
                    flag=1;
                    finished_one_component=1;
                    break;/*curve cannot be grown further*/
                }
        }
    }
    glutReshapeFunc(reshape);
    glutKeyboardFunc(kbd);
    glutDisplayFunc(display);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glColor3f(0.0, 0.0, 0.0);
    glLineWidth(3.0);
    glutMainLoop();
    exit(0);
    return 0;
}