/*
  Modified from An Introduction to OpenGL Programming,
  Ed Angel and Dave Shreiner, SIGGRAPH 2013
  Modified from 03_colorcube_rotate tutorial by Parag Chaudhuri, 2015
*/


#include "draft1.hpp"
#include <string.h>
#include <fstream>
#include <iostream>
using namespace std;

GLuint shaderProgram;
GLuint vbo, vao;

glm::mat4 rotation_matrix;
glm::mat4 ortho_matrix;
glm::mat4 modelview_matrix;
GLuint uModelViewMatrix;

//Skeleton data storage
vector<vector<int>> history,E,stage2junct;
vector<int> stage2nodes;
vector<vector<double>> V,U;
double xpos1,ypos1,xpos2,ypos2;
int fetchedpoints = 0;
//-----------------------------------------------------------------

//6 faces, 2 triangles/face, 3 vertices/triangle
const int num_vertices = 5000;

int tri_idx=0;
bool fetched_first=false;
int srcnodeindex=0,destnodeindex=0;
glm::vec4 v_positions[num_vertices];
glm::vec4 v_colors[num_vertices];

glm::vec4 avail_colors[4] = {
    glm::vec4(1.0f,0.0f,0.0f,1.0f),
    glm::vec4(0.0f,1.0f,0.0f,1.0f),
    glm::vec4(0.0f,0.0f,1.0f,1.0f),
    glm::vec4(1.0f,1.0f,0.0f,1.0f)
};

void loadvertexdata()
{
    FILE *fp;
    int i;
    fp = fopen("vertexdata.txt","r");
    fscanf(fp,"%d\n",&i);
    double x,y;
    while(i>0)
    {
        if(fscanf(fp,"%lf %lf\n",&x,&y)==2)
        {
            vector<double> temp;
            temp.push_back(x);
            temp.push_back(y);
            V.push_back(temp);
            i--;
        }
    }
    fscanf(fp,"%d\n",&i);
    while(i>0)
    {
        if(fscanf(fp,"%lf %lf\n",&x,&y)==2)
        {
            vector<double> temp;
            temp.push_back(x);
            temp.push_back(y);
            U.push_back(temp);
            i--;
        }
    }
    fscanf(fp,"%d\n",&i);
    int dumx,dumy;
    while(i>0)
    {
        if(fscanf(fp,"%d %d\n",&dumx,&dumy)==2)
        {
            vector<int> temp;
            temp.push_back(dumx);
            temp.push_back(dumy);
            E.push_back(temp);
            i--;
        }
    }
}

void loadSkeletondata()
{
    FILE *fp;
    int i,nodex;
    fp = fopen("skeldata.txt","r");
    fscanf(fp,"%d\n",&i);
    while(i>0)
    {
      fscanf(fp,"%d ",&nodex);
      stage2nodes.push_back(nodex);
      i--;
    }
    fscanf(fp,"\n");
    fscanf(fp,"%d\n",&i);
    while(i>0)
    {
      int hsize;
      fscanf(fp,"%d\n",&hsize);
      vector<int> temp;
      while(hsize>0)
      {
          fscanf(fp,"%d ",&nodex);
          temp.push_back(nodex);hsize--;
        }
        fscanf(fp,"\n");
        history.push_back(temp);
        i--;
    }
}

// generate 12 triangles: 36 vertices and 36 colors
void loadSkeleton(void)
{
  FILE *fp;
  int i;
  double x,y,z;
  // fp = fopen("stdskel.raw","r");
  fp = fopen("skeleton.raw","r");
  // fp = fopen("stand.raw","r");
  if(fp == NULL)
  {
     printf("Unable to open file for reading\n");
     exit(1);
  }
  i=0;
  while(fscanf(fp,"%lf %lf %lf\n",&x,&y,&z)==3)
  {
      cout<<x/120-1.0<<"  "<<-y/120 + 1.0<<endl;
      // v_positions[i] = glm::vec4(x/10,y/10,z/10,1.0f);
      v_positions[i]=glm::vec4(x/120 - 1.0, -y/120 + 1.0,z/120,1.0f);
      // v_positions[i]=glm::vec4(6*x/100, 6*y/100, 6*z/100, 1.0f);
      v_colors[i] = avail_colors[i%4];//glm::vec4(1.0f,1.0f,0.0f,1.0f);
      i++;
  }
  // num_vertices = i;
}

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && mode == 'M') // && mods == GLFW_MOD_SHIFT)
    {
        if(!fetched_first)
        {
            glfwGetCursorPos(window, &xpos1, &ypos1);
            GLdouble pos3D_x, pos3D_y, pos3D_z;
            GLdouble model_view[16];
glGetDoublev(GL_MODELVIEW_MATRIX, model_view);

GLdouble projection[16];
glGetDoublev(GL_PROJECTION_MATRIX, projection);

GLint viewport[4];
glGetIntegerv(GL_VIEWPORT, viewport);
gluUnProject(xpos1, ypos1, 0.01,model_view, projection, viewport,	&pos3D_x, &pos3D_y, &pos3D_z);
              cout<<pos3D_x<<" "<<pos3D_y<<" "<<pos3D_z<<endl;
            //Convert it to world coordinates
            //For now assuming that xpos1 and ypos1 are world coordinates
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U.at(stage2nodes.at(i)).at(0)/120 - 1.0;
                curry = U.at(stage2nodes.at(i)).at(1)/120 - 1.0;
                //pixel size is 10.. without convertion,
                if(xpos1 >= currx-5 && xpos1<= currx+5 && ypos1 >= curry-5 && ypos1<=curry-5)
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                //keep the -y in mind and make changes if needed
                {
                    srcnodeindex = stage2nodes.at(i);
                    fetched_first = !fetched_first; break;
                }
            }
        }
        else
        {
            glfwGetCursorPos(window, &xpos2, &ypos2);
            cout<<"x: "<<xpos2<<"  y: "<<ypos2<<endl;
            //convert it to world coordinates
            //Assuming the conversion done,
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U.at(stage2nodes.at(i)).at(0)/120 - 1.0;
                curry = U.at(stage2nodes.at(i)).at(1)/120 - 1.0;
                //pixel size is 10.. without convertion,
                if(xpos2 >= currx-5 && xpos2<= currx+5 && ypos2 >= curry-5 && ypos2<=curry-5 && i!=srcnodeindex)
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                //keep the -y in mind and make changes if needed
                {
                    destnodeindex = stage2nodes.at(i);
                    fetched_first = !fetched_first;
                    cout<<"fetched both the points. right click mouse to merge"<<endl;break;
                }
            }
        }
    }
    else if(button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS && mode == 'M' && !fetched_first && srcnodeindex != destnodeindex)
    {
        //Merging the srcnode and destnodeindex
        //Making vertex position as 0,0
        U.at(srcnodeindex).at(0) = 0.0; U.at(srcnodeindex).at(1) = 0.0;

        //Correcting vpositions and vcolors matrix data along with edge data.
        int positionindex = 0;
        //Removing the edges
        for(int i=0;i<E.size();i++)
        {
            if(!(E[i][0] == 0 && E[i][1] == 0))
            {
                if(E[i][0] == srcnodeindex && E[i][1] != destnodeindex)
                {
                  E[i][0] = destnodeindex;
                  double x,y;
                  x = U[destnodeindex][0]; y = U[destnodeindex][1];
                  v_positions[2*positionindex]=glm::vec4(x/120 - 1.0, -y/120 + 1.0,0.0f,1.0f);
                }
                else if(E[i][1] == srcnodeindex && E[i][0] != destnodeindex)
                {
                  E[i][1] = destnodeindex;
                  double x,y;
                  x = U[destnodeindex][0]; y = U[destnodeindex][1];
                  v_positions[2*positionindex + 1]=glm::vec4(x/120 - 1.0, -y/120 + 1.0,0.0f,1.0f);
                }
                else if((E[i][0] == srcnodeindex && E[i][1] == destnodeindex)|| (E[i][1] == srcnodeindex && E[i][0] == destnodeindex))
                {
                  E[i][1] = 0;  E[i][0] = 0;
                  for(int ki = 2*positionindex;ki<num_vertices;ki++)
                  {
                      v_positions[ki] = v_positions[ki+2];
                      v_colors[ki] = v_colors[ki+2];
                  }
                }
                positionindex++;
            }
        }

        //change history of destination node
        for(int i=0;i<history.at(srcnodeindex).size();i++)
        {
            if(find(history.at(destnodeindex).begin(),history.at(destnodeindex).end(),history.at(srcnodeindex).at(i)) == history.at(destnodeindex).end())
            {
                history.at(destnodeindex).push_back(history.at(srcnodeindex).at(i));
            }
        }

        //Remove from nodeslist;
        int indextoremove;
        indextoremove = find(stage2nodes.begin(),stage2nodes.end(),srcnodeindex)-stage2nodes.begin();
        stage2nodes.erase(stage2nodes.begin()+indextoremove);
        cout<<"Merging done"<<endl;
        srcnodeindex = 0; destnodeindex = 0;
    }
    else if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS && mode == 'A')
    {
        if(fetchedpoints == 0)
        {
            glfwGetCursorPos(window, &xpos1, &ypos1);
            //Convert it to world coordinates
            //For now assuming that xpos1 and ypos1 are world coordinates
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U.at(stage2nodes.at(i)).at(0)/120 - 1.0;
                curry = U.at(stage2nodes.at(i)).at(1)/120 - 1.0;
                //pixel size is 10.. without convertion,
                if(xpos1 >= currx-5 && xpos1<= currx+5 && ypos1 >= curry-5 && ypos1<=curry-5)
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                //keep the -y in mind and make changes if needed
                {
                    srcnodeindex = stage2nodes.at(i);
                    fetchedpoints = (fetchedpoints+1)%3; break;
                }
            }
        }
        else if(fetchedpoints == 1)
        {
            glfwGetCursorPos(window, &xpos2, &ypos2);
            cout<<"x: "<<xpos2<<"  y: "<<ypos2<<endl;
            //convert it to world coordinates
            //Assuming the conversion done,
            for(int i=0;i<stage2nodes.size();i++)
            {
                double currx,curry;
                currx = U.at(stage2nodes.at(i)).at(0)/120 - 1.0;
                curry = U.at(stage2nodes.at(i)).at(1)/120 - 1.0;
                //pixel size is 10.. without convertion,
                if(xpos2 >= currx-5 && xpos2<= currx+5 && ypos2 >= curry-5 && ypos2<=curry-5 && i!=srcnodeindex)
                //If the clicked position is somewhere in the 10*10 square of the node, we pick the node as src
                //keep the -y in mind and make changes if needed
                {
                    destnodeindex = stage2nodes.at(i);
                    fetchedpoints = (fetchedpoints+1)%3;
                    cout<<"fetched both the points. Select a position in between them"<<endl;break;
                }
            }
        }
        else
        {
            glfwGetCursorPos(window, &xpos1, &ypos1);
            //Convert it to world coordinates
            //For now assuming that xpos1 and ypos1 are world coordinates
            double xdist1,xdist2,ydist1,ydist2;
            xdist1 = xpos1 - U[srcnodeindex][0]; xdist1 = xdist1 > 0? xdist1 : -1 * xdist1;
            ydist1 = ypos1 - U[srcnodeindex][1]; ydist1 = ydist1 > 0? ydist1 : -1 * ydist1;
            xdist2 = xpos1 - U[destnodeindex][0]; xdist2 = xdist2 > 0? xdist2 : -1 * xdist2;
            ydist2 = ypos1 - U[destnodeindex][1]; ydist2 = ydist2 > 0? ydist2 : -1 * ydist2;

            double slopeofedge = (U[destnodeindex][1] - U[srcnodeindex][1]) / (U[destnodeindex][0] - U[srcnodeindex][0]);

            if(xdist1 > ydist1 && xdist2 > ydist2)
            //y coordinate of new point to be changed
            {
                ypos1 = (slopeofedge * (xpos1 - U[srcnodeindex][0])) + U[srcnodeindex][1];
            }
            else if(ydist1 > xdist1 && ydist2 > xdist2)
            //x coordinate to be changed
            {
                xpos1 = ((ypos1 - U[srcnodeindex][1]) / slopeofedge) + U[srcnodeindex][0];
            }

            vector<double> temp;
            temp.push_back(xpos1); temp.push_back(ypos1);

            //Insert in original vertices and in modified vertices
            V.push_back(temp);
            U.push_back(temp);

            //index of new node to be used in edges and dneighborhood, history
            int newindex = V.size() - 1;

            //Adding an entry in history.
            vector<int> dummy;
            dummy.push_back(newindex);
            history.push_back(dummy);

            int positionindex = 0;int ki;
            for(int i=0;i<E.size();i++)
            {
                if(!(E[i][0]==0 && E[i][1]))
                {
                    if((E[i][0] == srcnodeindex && E[i][1] == destnodeindex) || (E[i][0] == destnodeindex && E[i][1] == srcnodeindex))
                    {
                        E[i][0] = 0; E[i][1] = 0;
                        for(ki = 2*positionindex;ki<num_vertices;ki++)
                        {
                            v_positions[ki] = v_positions[ki+2];
                            v_colors[ki] = v_colors[ki+2];
                        }
                        break;
                    }
                    positionindex++;
                }
            }

            v_positions[ki] = glm::vec4((U[srcnodeindex][0])/120 - 1.0, -(U[srcnodeindex][1])/120 + 1.0,0.0f,1.0f);
            v_colors[ki] = avail_colors[ki%4]; ki++;
            v_positions[ki] = glm::vec4((U[newindex][0])/120 - 1.0, -(U[newindex][1])/120 + 1.0,0.0f,1.0f);
            v_colors[ki] = avail_colors[ki%4]; ki++;
            v_positions[ki] = glm::vec4((U[newindex][0])/120 - 1.0, -(U[newindex][1])/120 + 1.0,0.0f,1.0f);
            v_colors[ki] = avail_colors[ki%4]; ki++;
            v_positions[ki] = glm::vec4((U[destnodeindex][0])/120 - 1.0, -(U[destnodeindex][1])/120 + 1.0,0.0f,1.0f);
            v_colors[ki] = avail_colors[ki%4]; ki++;

            //Add new edges to edges list
            vector<int> temp1;
            temp1.push_back(srcnodeindex); temp1.push_back(newindex); E.push_back(temp1);
            temp1.push_back(newindex); temp1.push_back(destnodeindex); E.push_back(temp1);

            //Adjusting History of src and new node
            for(int i=0; i<history.at(srcnodeindex).size(); i++)
            {
                int nodeindexiter = history[srcnodeindex][i];
                //distance from node to src
                double dist1 = sqrt(pow((U[srcnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[srcnodeindex][1] - V[nodeindexiter][1]),2));
                double dist2 = sqrt(pow((U[newindex][0] - V[nodeindexiter][0]),2)+pow((U[newindex][1] - V[nodeindexiter][1]),2));
                if(dist1 > dist2) // Node is farther from original node than it is from new node
                {
                    history.at(newindex).push_back(nodeindexiter);
                    history.at(srcnodeindex).erase(history.at(srcnodeindex).begin()+i);
                }
            }

            //Adjusting History of src and new node
            for(int i=0; i<history.at(destnodeindex).size(); i++)
            {
                int nodeindexiter = history[destnodeindex][i];
                //distance from node to src
                double dist1 = sqrt(pow((U[destnodeindex][0] - V[nodeindexiter][0]),2)+pow((U[destnodeindex][1] - V[nodeindexiter][1]),2));
                double dist2 = sqrt(pow((U[newindex][0] - V[nodeindexiter][0]),2)+pow((U[newindex][1] - V[nodeindexiter][1]),2));
                if(dist1 > dist2) // Node is farther from original node than it is from new node
                {
                    history.at(newindex).push_back(nodeindexiter);
                    history.at(destnodeindex).erase(history.at(destnodeindex).begin()+i);
                }
            }
            fetchedpoints = (fetchedpoints+1)%3;
        }
    }
}

void initBuffersGL(void)
{
  //cg::loadSkeleton();
  // loadvertexdata(); //to debug.
  loadSkeletondata();
  loadSkeleton();

  //Ask GL for a Vertex Attribute Object (vao)
  glGenVertexArrays (1, &vao);
  //Set it as the current array to be used by binding it
  glBindVertexArray (vao);

  //Ask GL for a Vertex Buffer Object (vbo)
  glGenBuffers (1, &vbo);
  //Set it as the current buffer to be used by binding it
  glBindBuffer (GL_ARRAY_BUFFER, vbo);
  //Copy the points into the current buffer
  glBufferData (GL_ARRAY_BUFFER, sizeof (v_positions) + sizeof(v_colors), NULL, GL_DYNAMIC_DRAW);
  glBufferSubData( GL_ARRAY_BUFFER, 0, sizeof(v_positions), v_positions );
  glBufferSubData( GL_ARRAY_BUFFER, sizeof(v_positions), sizeof(v_colors), v_colors );

  // Load shaders and use the resulting shader program
  std::string vertex_shader_file("03_vshader.glsl");
  std::string fragment_shader_file("03_fshader.glsl");

  std::vector<GLuint> shaderList;
  shaderList.push_back(csX75::LoadShaderGL(GL_VERTEX_SHADER, vertex_shader_file));
  shaderList.push_back(csX75::LoadShaderGL(GL_FRAGMENT_SHADER, fragment_shader_file));

  shaderProgram = csX75::CreateProgramGL(shaderList);
  glUseProgram( shaderProgram );

  // set up vertex arrays
  GLuint vPosition = glGetAttribLocation( shaderProgram, "vPosition" );
  glEnableVertexAttribArray( vPosition );
  glVertexAttribPointer( vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0) );

  GLuint vColor = glGetAttribLocation( shaderProgram, "vColor" );
  glEnableVertexAttribArray( vColor );
  glVertexAttribPointer( vColor, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(sizeof(v_positions)) );

  uModelViewMatrix = glGetUniformLocation( shaderProgram, "uModelViewMatrix");
}

void renderGL(void)
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  rotation_matrix = glm::rotate(glm::mat4(1.0f), xrot, glm::vec3(1.0f,0.0f,0.0f));
  rotation_matrix = glm::rotate(rotation_matrix, yrot, glm::vec3(0.0f,1.0f,0.0f));
  rotation_matrix = glm::rotate(rotation_matrix, zrot, glm::vec3(0.0f,0.0f,1.0f));
  ortho_matrix = glm::ortho(-2.0, 2.0, -2.0, 2.0, -2.0, 2.0);

  modelview_matrix = ortho_matrix * rotation_matrix;

  glUniformMatrix4fv(uModelViewMatrix, 1, GL_FALSE, glm::value_ptr(modelview_matrix));
  // Draw
  glDrawArrays(GL_LINES, 0, num_vertices);
  glDrawArrays(GL_POINTS, 0, num_vertices);

}

int main(int argc, char** argv)
{
  //! The pointer to the GLFW window
  GLFWwindow* window;

  //! Setting up the GLFW Error callback
  glfwSetErrorCallback(csX75::error_callback);

  //! Initialize GLFW
  if (!glfwInit())
    return -1;

  //We want OpenGL 4.0
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  //This is for MacOSX - can be omitted otherwise
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  //We don't want the old OpenGL
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  //! Create a windowed mode window and its OpenGL context
  window = glfwCreateWindow(900, 900, "Assignment 1", NULL, NULL);
  if (!window)
  {
    glfwTerminate();
    return -1;
  }

  //! Make the window's context current
  glfwMakeContextCurrent(window);

  //Initialize GLEW
  //Turn this on to get Shader based OpenGL
  glewExperimental = GL_TRUE;
  GLenum err = glewInit();
  if (GLEW_OK != err)
    {
      //Problem: glewInit failed, something is seriously wrong.
      std::cerr<<"GLEW Init Failed : %s"<<std::endl;
    }

  //Print and see what context got enabled
  std::cout<<"Vendor: "<<glGetString (GL_VENDOR)<<std::endl;
  std::cout<<"Renderer: "<<glGetString (GL_RENDERER)<<std::endl;
  std::cout<<"Version: "<<glGetString (GL_VERSION)<<std::endl;
  std::cout<<"GLSL Version: "<<glGetString (GL_SHADING_LANGUAGE_VERSION)<<std::endl;

  //Keyboard Callback
  glfwSetKeyCallback(window, csX75::key_callback);
  //Framebuffer resize callback
  glfwSetFramebufferSizeCallback(window, csX75::framebuffer_size_callback);
  glfwSetMouseButtonCallback(window, mouse_button_callback);

  // Ensure we can capture the escape key being pressed below
  glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_TRUE);

  //Initialize GL state
  csX75::initGL();
  initBuffersGL();

  // Loop until the user closes the window
  while (glfwWindowShouldClose(window) == 0)
    {

      // Render here
      renderGL();

      // Swap front and back buffers
      glfwSwapBuffers(window);

      // Poll for and process events
      glfwPollEvents();
    }

  glfwTerminate();
  return 0;
}

//-------------------------------------------------------------------------
