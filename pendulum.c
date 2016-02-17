#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <math.h>
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xlocale.h>

#define N 10000
#define L1 1.4
#define L2 1.2
#define M1 1.2
#define M2 0.4
#define THETA1 3.5
#define THETA2 1.0
#define GRAV 9.8

typedef struct __point {
  double x;
  double y;
} Point;

typedef struct __pendulum {
  double theta1;
  double theta2;
  double omega1;
  double omega2;
  Point p1_p; // 1つ目の振り子のデカルト座標系における座標
  Point p2_p; // 2つ目　〃
} Pendulum;

const int BORDER = 2;
const int WIDTH  = 640;
const int HEIGHT = 480;
unsigned long black, white, red, blue;
double l1;   // 1つ目の振り子の長さ
double l2;   // 2つ目 〃
double m1;   // 1つ目の振り子の質量
double m2;   // 2つ目 〃
const double g = GRAV;  // 重力加速度
Pendulum p;


double f1(double theta1, double omega1, double theta2, double omega2) {
  // 角加速度の計算
  return (-(m1+m2)*g*sin(theta1)/l1 + m2*g*cos(theta1-theta2)*sin(theta2)/l1
          -m2*cos(theta1-theta2)*sin(theta1-theta2)*omega1*omega1 - m2*sin(theta1-theta2)*omega2*omega2*(l2/l1))
         / (m1+m2*sin(theta1-theta2)*sin(theta1-theta2));
}

double f2(double theta1, double omega1, double theta2, double omega2) {
  // 角加速度の計算
  return ((m1+m2)*g*cos(theta1)*sin(theta1-theta2)/l2 + (m1+m2)*sin(theta1-theta2)*omega1*omega1*(l1/l2)
          + m2*cos(theta1-theta2)*sin(theta1-theta2)*omega2*omega2)
         / (m1+m2*sin(theta1-theta2)*sin(theta1-theta2));
}

double g1(double omega1) {
  return omega1;  // 角速度
}

double g2(double omega2) {
  return omega2;  // 角速度
}


// 4次のルンゲクッタ法
void Runge_Kutta(void) {
  int step;
  double k1[4], k2[4], k3[4], k4[4];
  double dt = 0.000001;

  for (step = 0; step < N; step++) {
    k1[0] = dt * g1(p.omega1);
    k1[1] = dt * f1(p.theta1, p.omega1, p.theta2, p.omega2);
    k1[2] = dt * g2(p.omega2);
    k1[3] = dt * f2(p.theta1, p.omega1, p.theta2, p.omega2);

    k2[0] = dt * g1(p.omega1+k1[1]/2.0);
    k2[1] = dt * f1(p.theta1+k1[0]/2.0, p.omega1+k1[1]/2.0, p.theta2+k1[2]/2.0, p.omega2+k1[3]/2.0);
    k2[2] = dt * g2(p.omega2+k1[3]/2.0);
    k2[3] = dt * f2(p.theta1+k1[0]/2.0, p.omega1+k1[1]/2.0, p.theta2+k1[2]/2.0, p.omega2+k1[3]/2.0);

    k3[0] = dt * g1(p.omega1+k2[1]/2.0);
    k3[1] = dt * f1(p.theta1+k2[0]/2.0, p.omega1+k2[1]/2.0, p.theta2+k2[2]/2.0, p.omega2+k2[3]/2.0);
    k3[2] = dt * g2(p.omega2+k2[3]/2.0);
    k3[3] = dt * f2(p.theta1+k2[0]/2.0, p.omega1+k2[1]/2.0, p.theta2+k2[2]/2.0, p.omega2+k2[3]/2.0);

    k4[0] = dt * g1(p.omega1+k3[1]);
    k4[1] = dt * f1(p.theta1+k3[0], p.omega1+k3[1], p.theta2+k3[2], p.omega2+k3[3]);
    k4[2] = dt * g2(p.omega2+k3[3]);
    k4[3] = dt * f2(p.theta1+k3[0], p.omega1+k3[1], p.theta2+k3[2], p.omega2+k3[3]);

    p.theta1 += (k1[0] + 2*k2[0] + 2*k3[0] + k4[0]) / 6.0;
    p.omega1 += (k1[1] + 2*k2[1] + 2*k3[1] + k4[1]) / 6.0;
    p.theta2 += (k1[2] + 2*k2[2] + 2*k3[2] + k4[2]) / 6.0;
    p.omega2 += (k1[3] + 2*k2[3] + 2*k3[3] + k4[3]) / 6.0;

    p.p2_p.x = (l1 * sin(p.theta1) + l2 * sin(p.theta2));
    p.p2_p.y = (l1 * cos(p.theta1) + l2 * cos(p.theta2));
    p.p1_p.x = l1 * sin(p.theta1);
    p.p1_p.y = l1 * cos(p.theta1);
  }
}

// x軸・y軸の描写
void drawAxis(Display *dpy, Window w, GC gc) {
    XSetForeground(dpy, gc, black);
    XDrawLine(dpy, w, gc, WIDTH/2, 10, WIDTH/2, HEIGHT-10);
    XDrawLine(dpy, w, gc, 10, HEIGHT/2, WIDTH-10, HEIGHT/2);
    XDrawString(dpy, w, gc, WIDTH/2-10, HEIGHT/2+16, "O", 1);
}


int main(int argc, char const *argv[]) {
  Display *dpy;
  int screen;
  Window root, w;
  GC gc;
  XEvent e;
  FILE *fp1, *fp2;
  Colormap cmap;
  XColor color, exact;

  if (argc == 1) {
    l1 = L1;
    l2 = L2;
    m1 = M1;
    m2 = M2;
    p.theta1 = THETA1;
    p.theta2 = THETA2;
  } else if (argc == 7) {
    l1 = atof(argv[1]);
    l2 = atof(argv[2]);
    m1 = atof(argv[3]);
    m2 = atof(argv[4]);
    p.theta1 = atof(argv[5]);
    p.theta2 = atof(argv[6]);
  } else {
    fprintf(stderr, "Error invalid parameter.\n");
    fprintf(stderr, "Usage: \"$ ./pendulum l1 l2 m1 m2 theta1 theta2\"\n       \"$ ./pendulum\" (Initial Parameter)\n");
    exit(1);
  }

  p.omega1 = 0.0;
  p.omega2 = 0.0;
  p.p1_p.x = l1 * sin(p.theta1);
  p.p1_p.y = l1 * cos(p.theta1);
  p.p2_p.x = (l1 * sin(p.theta1) + l2 * sin(p.theta2));
  p.p2_p.y = (l1 * cos(p.theta1) + l2 * cos(p.theta2));

  setlocale(LC_ALL, "");
  dpy = XOpenDisplay("");

  root = DefaultRootWindow(dpy);
  screen = DefaultScreen(dpy);

  cmap = DefaultColormap(dpy, screen);

  white = WhitePixel(dpy, screen);
  black = BlackPixel(dpy, screen);
  XAllocNamedColor( dpy, cmap, "red", &color, &exact);
  red = color.pixel;
  XAllocNamedColor( dpy, cmap, "blue", &color, &exact);
  blue = color.pixel;

  w = XCreateSimpleWindow(dpy, root, 100, 100, WIDTH, HEIGHT, BORDER, black, white);
  gc = XCreateGC(dpy, w, 0, NULL);

  XSelectInput(dpy, w, ButtonPressMask);
  XMapWindow(dpy, w);

  fp1 = fopen("pendulum_m1.out", "w");
  fp2 = fopen("pendulum_m2.out", "w");

  XSetForeground(dpy, gc, black);

  while(1) {
    XSetForeground( dpy, gc, white);
    XDrawLine(dpy, w, gc, p.p1_p.x, p.p1_p.y, WIDTH/2, HEIGHT/2);
    XDrawLine(dpy, w, gc, p.p2_p.x, p.p2_p.y, p.p1_p.x, p.p1_p.y);
    XFillArc(dpy, w, gc, p.p2_p.x-5, p.p2_p.y-5, 10, 10, 0, 360*64);
    XFillArc(dpy, w, gc, p.p1_p.x-6, p.p1_p.y-6, 12, 12, 0, 360*64);
    XSetForeground(dpy, gc, black);

    drawAxis(dpy, w, gc);

    Runge_Kutta();

    p.p1_p.x *= 100, p.p1_p.y *= 100;
    p.p1_p.x += WIDTH/2, p.p1_p.y += HEIGHT/2;
    XDrawLine(dpy, w, gc, p.p1_p.x, p.p1_p.y, WIDTH/2, HEIGHT/2);
    p.p2_p.x *= 100, p.p2_p.y *= 100;
    p.p2_p.x += WIDTH/2, p.p2_p.y += HEIGHT/2;
    XDrawLine(dpy, w, gc, p.p2_p.x, p.p2_p.y, p.p1_p.x, p.p1_p.y);

    XSetForeground( dpy, gc, red);
    XFillArc(dpy, w, gc, p.p2_p.x-5, p.p2_p.y-5, 10, 10, 0, 360*64);
    XSetForeground( dpy, gc, blue);
    XFillArc(dpy, w, gc, p.p1_p.x-6, p.p1_p.y-6, 12, 12, 0, 360*64);

    fprintf(fp1, "%lf %lf\n", -(p.p1_p.x-WIDTH/2)/100, -(p.p1_p.y-HEIGHT/2)/100);
    fprintf(fp2, "%lf %lf\n", -(p.p2_p.x-WIDTH/2)/100, -(p.p2_p.y-HEIGHT/2)/100);

    if (!XPending(dpy)) {
        usleep(10);
        continue;
    }
  }

  fclose(fp1);
  fclose(fp2);
  return 0;
}
