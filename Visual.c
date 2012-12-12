#include <math.h>
#include <stdio.h>
#include <GL/glut.h>


void draw(void);
void drawElement(float x1, float y1, float x2, float y2, float node1Fac, float node2Fac);
void drawNode(float x, float y);
void drawResultValues(int mode);
void drawSystem(void);

GLuint Window = 0;
float scale, center_x, center_y, max_N, max_V, max_M;
float **elements;
int nodes, size, mode; //Mode = 8 displays normal forces, Mode = 9 displays shear forces, Mode = 10 displays moments.

void readFile(void) {
  FILE *file;
  int i;
  char line[500];
  float length, angle, node1, node2, garbage, u1, v1, w1, u2, v2, w2;
  
  file = fopen ("input1.dat", "rt");
  fgets(line, 500, file);                            //Reads 'Number of elements' line and do nothing.
  fscanf(file,"%d",&size);                           //Reads number of elements.
  
  fgets(line, 500, file);                            //Reads second line and do nothing.
  fgets(line, 500, file);                            //Reads 'Number of nodes' line and do nothing.
  fscanf(file, "%d", &nodes);                        //Reads number of nodes.
  fgets(line, 500, file);                            //Reads 'Element properties' line and do nothing.
  fgets(line,500,file);                              //Reads first element row.
  elements = (float**)malloc(size * sizeof(float*)); //Allocate matrix of size size*14. 
  for(i=0;i<size; i++) {
    elements[i] = (float*)malloc(14 * sizeof(float));
  }
  for(i=0; i<size; i++){
    fscanf(file,"%f %f %f %f %f %f %f", &garbage, &garbage, &garbage, &length, &angle, &node1, &node2);
    elements[i][0] = length;
    elements[i][1] = angle;
    elements[i][2] = node1;
    elements[i][3] = node2;
    fgets(line,500,file);
  }
  fclose(file);

  file = fopen("result.dat", "rt");
  fgets(line, 500, file); //Reads 'RESULT FILE:' line and do nothing.
  fgets(line, 500, file); //Reads 'System displacements' line and do nothing.
  fgets(line, 500, file); //Reads displacements line and do nothing.
  fgets(line, 500, file); //Reads 'Nodal forces of each element' line and do nothing.
  for(i=0; i<size; i++){
    fscanf(file, "%f %f %f %f %f %f", &u1, &v1, &w1, &u2, &v2, &w2);
    elements[i][8] = u1;
    elements[i][9] = v1;
    elements[i][10] = w1;
    elements[i][11] = u2;
    elements[i][12] = v2;
    elements[i][13] = w2;
    fgets(line,500,file);
  }
  fclose(file); 
}

/*
 * Constructs a physical coordinate for every node.
 */
void constructNodeCoordinates(void) {
  float x_2, y_2, u, v, length, angle;
  int i,j, found;
  float visited[nodes][3];

  for(i=0;i<size;i++) {
    length = elements[i][0];
    angle = elements[i][1];
    u = length * cos(angle);
    v = length * sin(angle);
    for(j=0;j<nodes;j++) {
      if(elements[i][2] == visited[j][0]) {
	elements[i][4] = visited[j][1];
	elements[i][5] = visited[j][2];
      }
    }
    x_2 = elements[i][4] + u;
    y_2 = elements[i][5] + v;
    elements[i][6] = x_2;
    elements[i][7] = y_2;
    visited[(int)elements[i][3]][0] = elements[i][3];
    visited[(int)elements[i][3]][1] = x_2;
    visited[(int)elements[i][3]][2] = y_2;
  }
}

void scaleSystem(void) {
  float max_x, max_y, min_x, min_y, scale;
  int i;

  for(i=0;i<size;i++) {
    if(elements[i][4] > max_x)
      max_x = elements[i][4];
    if(elements[i][5] > max_y)
      max_y = elements[i][5];
    if(elements[i][4] < min_x)
      min_x = elements[i][4];
    if(elements[i][5] < min_y)
      min_y = elements[i][5];
    if(elements[i][6] > max_x)
      max_x = elements[i][6];
    if(elements[i][7] > max_y)
      max_y = elements[i][7];
    if(elements[i][6] < min_x)
      min_x = elements[i][6];
    if(elements[i][7] < min_y)
      min_y = elements[i][7];
  }

  if((max_x - min_x) > (max_y - min_y))
    scale = max_x - min_x;
  else
    scale = max_y - min_y;

  center_x = ((max_x + min_x) / 2) / scale;
  center_y = ((max_y + min_y) / 2) / scale;

  for(i=0;i<size;i++) {
    elements[i][4] = elements[i][4] / scale;
    elements[i][5] = elements[i][5] / scale;
    elements[i][6] = elements[i][6] / scale;
    elements[i][7] = elements[i][7] / scale;
  }
}

/*
 * This procedure is not in use.
 */
void accumulateNodeForces(void) {
  int i,j;
  float accumulated;

  for(i=0;i<size;i++) {
    for(j=0;j<size;j++) {
      if(elements[i][4] == elements[i][6] && elements[j][5] == elements[j][7]) { //Accumulate normal forces.
	accumulated = elements[i][8] + elements[j][11];
	elements[i][8] = accumulated;
	elements[j][11] = accumulated;
      }
      else if(elements[i][4] == elements[i][6] && elements[j][5] == elements[j][7]) { //Accumulate shear forces.
	accumulated = elements[i][9] + elements[j][12];
	elements[i][9] = accumulated;
	elements[j][12] = accumulated;
      }
      else if(elements[i][4] == elements[i][6] && elements[j][5] == elements[j][7]) { //Accumulate moment forces.
	accumulated = elements[i][10] + elements[j][12];
	elements[i][10] = accumulated;
	elements[j][12] = accumulated;
      }
    }
  }
}

/*
 * Finds maximum values of forces and scale other forces after the maximum.
 */
void constructNodeColor(void) {
  int i;
  float maxNormal, maxShear, maxMoment;

  for(i=0;i<size;i++) {
    if(fabsf(elements[i][8]) > maxNormal)
      maxNormal = fabsf(elements[i][8]);
    if(fabsf(elements[i][11]) > maxNormal)
      maxNormal = fabsf(elements[i][11]);
    if(fabsf(elements[i][9]) > maxShear) 
      maxShear = fabsf(elements[i][9]);
    if(fabsf(elements[i][12]) > maxShear) 
      maxShear = fabsf(elements[i][12]);
    if(fabsf(elements[i][10]) > maxMoment) 
      maxMoment = fabsf(elements[i][10]);
    if(fabsf(elements[i][13]) > maxMoment)
      maxMoment = fabsf(elements[i][13]);
  }

  max_N = maxNormal;
  max_V = maxShear;
  max_M = maxMoment;
  
  for(i=0;i<size;i++) {
    if(maxNormal != 0) {
      elements[i][8] = elements[i][8] / maxNormal;
      elements[i][11] = elements[i][11] / maxNormal;
    }
    if(maxShear != 0) {
      elements[i][9] = elements[i][9] / maxShear;
      elements[i][12] = elements[i][12] / maxShear;
    }
    if(maxMoment != 0) {
      elements[i][10] = elements[i][10] / maxMoment;
      elements[i][13] = elements[i][13] / maxMoment;
    }
  }
}

void draw(void) {
  
  glClearColor(1,1,1,0);
  glClear(GL_COLOR_BUFFER_BIT);
  glLoadIdentity();
  drawResultValues(mode);
  glTranslatef(-center_x,-center_y,0);
  drawSystem();
  glutSwapBuffers();
}

void drawSystem (void) {
  int i;
  for(i=0;i<size;i++) {
    drawElement(elements[i][4], elements[i][5], elements[i][6], elements[i][7], elements[i][mode], elements[i][mode+3]);
  }
}

void drawElement (float x1, float y1, float x2, float y2, float node1Fac, float node2Fac) {
  
  float red1 = 0.7 + 0.3 * fabsf(node1Fac);
  float red2 = 0.7 + 0.3 * fabsf(node2Fac);
  float other1 = 0.7 - 0.3 * fabsf(node1Fac);
  float other2 = 0.7 - 0.3 * fabsf(node2Fac);
  
  glLineWidth(30);
  glBegin(GL_LINES);
  glColor3f(red1,other1,other1); //Gray colored element
  glVertex2f(x1,y1);
  glColor3f(red2,other2,other2);
  glVertex2f(x2, y2);
  glEnd();
}

void drawNode (float x, float y) {
  float angle,x2,y2;
  float radius = 0.015;
  glBegin(GL_TRIANGLE_FAN);
  glColor3f(0,0,0);
  glVertex2f(x,y);

  for (angle=1.0f; angle<361.0f; angle+=0.2) {
    x2 = x+sin(angle)*radius;
    y2 = y+cos(angle)*radius;
    glVertex2f(x2,y2);
  }
  glEnd();
}

void drawResultValues(int mode) {
  char value[40];
  glColor3f(0,0,0);
  glRasterPos3f(-0.8,0.8,0);
  switch(mode) {
  case 8:
    sprintf(value,"%s %d %s","Max. normal force:", (int)max_N, "N");
    glutBitmapString(GLUT_BITMAP_HELVETICA_18, value);
    break;
  case 9:
    sprintf(value,"%s %d %s","Max. shear force:", (int)max_V, "N");
    glutBitmapString(GLUT_BITMAP_HELVETICA_18, value);
   break;
  case 10:
    sprintf(value,"%s %d %s","Max. moment:", (int)max_M, "Nmm");
    glutBitmapString(GLUT_BITMAP_HELVETICA_18, value);
    break;
  }
}

void key(unsigned char k, int x, int y) {
  (void) x;
  (void) y;
  switch (k) {
  case 27:
    glutDestroyWindow(Window);
    exit(0);
    break;
  case '1':
    mode = 8;
    drawSystem();
    glutPostRedisplay();
    break;
  case '2':
    mode = 9;
    drawSystem();
    glutPostRedisplay();
    break;
  case '3':
    mode = 10;
    drawSystem();
    glutPostRedisplay();
    break;
  }
}

int main(int argc, char **argv) {

  readFile();
  constructNodeCoordinates();
  scaleSystem();
  constructNodeColor();
  mode = 8;

  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
  glutInitWindowPosition(1000,25);
  glutInitWindowSize(900,700);
  glutCreateWindow("FEM simulation result");
  glutKeyboardFunc(key);
  glutDisplayFunc(draw);
  glutMainLoop();

  return 0;
}
