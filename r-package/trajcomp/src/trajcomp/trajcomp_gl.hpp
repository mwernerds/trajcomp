#ifndef TRAJCOMP_GL_INC
#define TRAJCOMP_GL_INC

/*OpenGL utility libraries for trajcomp*/
/*Note: these do not work on Apple Mac, though I dont care */
#include <GL/glu.h>
#include <GL/gl.h>

/*JET COLOR PALETTE*/
void glJet(double v, double vmin, double vmax, double alpha=1)
{
	double r=1, g =1, b = 1.0;
   double dv;
   if (v < 0)
   {
	  r = g = b = 1;
	  alpha = 0.1;
   }else{

   if (v < vmin)
      v = vmin;
   if (v > vmax)
      v = vmax;
   dv = vmax - vmin;

   if (v < (vmin + 0.25 * dv)) {
      r = 0;
      g = 4 * (v - vmin) / dv;
   } else if (v < (vmin + 0.5 * dv)) {
      r = 0;
      b = 1 + 4 * (vmin + 0.25 * dv - v) / dv;
   } else if (v < (vmin + 0.75 * dv)) {
      r = 4 * (v - vmin - 0.5 * dv) / dv;
      b = 0;
   } else {
      g = 1 + 4 * (vmin + 0.75 * dv - v) / dv;
      b = 0;
   }
	}
   glColor4f(r,g,b,alpha);
}


/*A simple graph scaling on bounding box*/
template<typename value_type>
void glRenderGraph(size_t x, size_t y, size_t width, size_t height,
				std::vector<std::pair<value_type,value_type>> &graph)
{ 
		
		glViewport(x,y,width,height);
		// 1.) Fill background
		glMatrixMode( GL_PROJECTION );
		glLoadIdentity();
		glOrtho(0,1,0,1,-10,10); 
		
		glMatrixMode(GL_MODELVIEW);
		glColor3f(1,1,0.8);
		glBegin(GL_QUADS);
			glVertex2f(0,0);
			glVertex2f(0,1);
			glVertex2f(1,1);
			glVertex2f(1,0);
		glEnd();
		//2. if graph is there, render it
		if (graph.size() != 0)
		{
			glMatrixMode( GL_PROJECTION );
			glLoadIdentity();
			double ymin = std::numeric_limits<double>::infinity();
			double ymax = 0;
			for (size_t i=0; i < graph.size(); i++)
			{
				if (graph[i].second > ymax) ymax = graph[i].second;
				if (graph[i].second < ymin) ymin = graph[i].second;
			}
			double offset = (ymax-ymin) / 10;
			
			glOrtho(graph[0].first,graph[graph.size()-1].first,ymin-offset,ymax+offset,-10,10); // y x rendering
			glMatrixMode(GL_MODELVIEW);
			// Render X-Axis
			/*glColor3f(0,1,0);
			glBegin(GL_LINES);
			glVertex2f(graph[0].first,ymin);
			glVertex2f(graph[graph.size()-1].first,ymin);
			glEnd();*/
			
			glColor3f(1,0,0);
			glBegin(GL_LINE_STRIP);
			for (size_t i=0; i < graph.size(); i++)
			{
				glVertex2f(graph[i].first,graph[i].second);
			}
			
			
			glEnd();
		}
};

#endif
