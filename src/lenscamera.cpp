#include "lenscamera.h"

#include "image.h"

using namespace std;

namespace CGL {


/****** Helpers ******/
  

// Extract the R, G, or B channel out of an RGBA color stored in a single 32bit integer
static uint32_t red_channel(uint32_t color) {
    return (255 & (color >> (0)));
}

static uint32_t green_channel(uint32_t color) {
    return (255 & (color >> (8)));
}

static uint32_t blue_channel(uint32_t color) {
    return (255 & (color >> (16)));
}

// Convert from millimeters to meters
static const double scale = .001;







/****** LensElement functions ******/


bool LensElement::pass_through(Ray &r, double &prev_ior) const {
  // Part 1 Task 1: Implement this. It takes r and passes it through this lens element.
  Vector3D hit_p = Vector3D(0,0,0);
  if (intersect(r, &hit_p)) {
    if (pow(hit_p.x,2)+pow(hit_p.y,2)<=pow(aperture/2.0,2)) {
      if (radius == 0) {
        prev_ior = ior;
        return true;
      } 
      if (refract(r, hit_p, prev_ior)) {
        prev_ior = ior;
        return true;
      }
    }
  }
  return false;
}

bool LensElement::test(const Ray& r, double& t1, double& t2) const {

  // TODO Part 1, task 4:
  // Implement ray - sphere intersection test.
  // Return true if there are intersections and writing the
  // smaller of the two intersection times in t1 and the larger in t2.
  double a = dot(r.d, r.d);
  Vector3D cen = Vector3D(0,0,center);
  double b = 2*dot(r.o-cen, r.d);
  double c = dot(r.o-cen, r.o-cen) - pow(radius,2);
  if (pow(b,2) - 4*a*c < 0) return false;
  if (a == 0) return false;
  double delta = pow(b,2) - 4*a*c;
  if (delta < 0) {
    delta = 0;
  }
  t1 = (-b-sqrt(delta))/2/a;
  t2 = (-b+sqrt(delta))/2/a;
  return true;
}

bool LensElement::intersect(const Ray &r, Vector3D *hit_p) const {
  // Part 1 Task 1: Implement this. It intersects r with this spherical lens elemnent 
  // (or aperture diaphragm). You'll want to reuse some sphere intersect code.
  if (radius == 0.f) { //is aperture element
    if (r.d.z != 0.f) {
      double t = (0.f-r.o.z)/r.d.z;
      *hit_p= r.o+t*r.d;  
      return true;
    }
    return false;
  }
  // if normal lens element
  double t1, t2, t;
  if (test(r, t1, t2)) {
      Vector3D hit_p1 = r.o+t1*r.d;
      Vector3D hit_p2 = r.o+t2*r.d;
      bool same_side1 = ((center-hit_p1.z) * radius >= 0);
      bool same_side2 = ((center-hit_p2.z) * radius >= 0);
      if ((t1>=0) && (t1 >= r.min_t) && (t1 <= r.max_t) && same_side1) {
        t = t1; 
      } else if ((t2>=0) && (t2 >= r.min_t) && (t2 <= r.max_t) && same_side2) {
        t = t2; 
      } else {
        return false;
      } 
      *hit_p= r.o+t*r.d;  
      return true;
  }
  return false;
}
bool LensElement::refract(Ray& r, const Vector3D& hit_p, const double& prev_ior) const {
  // Part 1 Task 1: Implement this. It refracts the Ray r with this lens element or 
  // does nothing at the aperture element.
  // You'll want to consult your refract function from the previous assignment.
  //  \  i /| 
  //     \  | n
  //       \|  n1=prev_ior
  //  ----------
  //        |\ n2=ior
  //        | \t
  // Vector form of shell's law:
  // http://graphics.stanford.edu/courses/cs148-10-summer/docs/2006--degreve--reflection_refraction.pdf
  Vector3D n = (hit_p-Vector3D(0,0,center)).unit();
  if (r.d.z * radius < 0) {
      n = (-1) * n;
  }
  Vector3D i = r.d.unit();
  double n1,n2;
  if (r.d.z < 0) {
    n1 = prev_ior;
    n2 = ior;
  } else {
    n2 = prev_ior;
    n1 = ior;
  }
  double cos_theta_i = -dot(n,i);
  double sin2theta_i = 1-pow(cos_theta_i,2);
  if (n1*sqrt(sin2theta_i) <= n2) {
    double sin2theta_t = pow(n1/n2,2)*sin2theta_i;
    Vector3D t = (n1/n2)*i+
                 (n1/n2*cos_theta_i-
                  sqrt(1-sin2theta_t))*n;

    r.d = t.unit();
    r.o = hit_p;
    return true;
  }
  return false;
}






/****** Lens functions ******/



void Lens::parse_lens_file(std::string filename) {

  ifstream infile(filename);
  string line;
  double z_coord = 0;
  double z_ap;
  vector<LensElement> backwards;
  elts.clear();
  bool first = true;
  while (getline(infile, line)) {
    if (first) {
      cout << "[Lens] Loading lens file " << line << endl;
      first = false;
    }
    if (line[0] == '#')
      continue;
    stringstream ss(line);
    LensElement lens;
    double offset;
    ss >> lens.radius >> offset >> lens.ior >> lens.aperture;
    lens.center = z_coord;
    if (!lens.radius) {
      z_ap = z_coord;
    }
    z_coord += offset;
    backwards.push_back(lens);
  }
  for (int i = backwards.size() - 1; i >= 0; --i) {
    LensElement l = backwards[i];
    l.center = (l.center - z_ap) + l.radius;
    if (i) l.ior = backwards[i-1].ior;
    else l.ior = 1;
    if (!l.ior) l.ior = 1;
    elts.push_back(l);
    if (!l.radius)
      ap_i = elts.size()-1;
    // cout << "Lens element edge first " << (l.center - l.radius) << " " 
    //   << l.radius << " " << l.center << " " << l.ior << " " << l.aperture << endl;
  }
  double c = elts.front().center, r = elts.front().radius, a = elts.front().aperture * .5;
  back_elt = c - (r>0?1:-1) * sqrt(r*r-a*a);
  ap_radius = ap_original = elts[ap_i].aperture;

  // Get infinity and close focus depths, also get focal length.
  set_focus_params();
  // Focus at infinity to start.
  sensor_depth = infinity_focus;
       
}

void Lens::set_focus_params() {

  // Part 1 Task 2: Implement this. 
  // After this function is called, the three variables
  // infinity_focus, near_focus, and focal_length
  // should be set correctly.
  /* ==== infinit_focus ========= */
  double epsilon = ap_radius/150.0;
  LensElement front_lens = elts.back();
  double front_lens_z = front_lens.center-front_lens.radius;
  double t_inf, t_near, t_focal;
  Ray r = Ray(Vector3D(epsilon, 0, front_lens_z-10), Vector3D(0,0,1));
  vector<Vector3D> trace;
  trace.push_back(r.o);
  if (trace_backwards(r, &trace)) {
      if ((r.o.x !=0) && (r.d.x == 0)) {
         t_inf = INF_D;
         infinity_focus = INF_D;
      } else {
         t_inf = -r.o.x/r.d.x;
         infinity_focus = r.o.z + t_inf*r.d.z;
      }
  } else {
      infinity_focus = -42;
  }
  /* ==== focal_focus ========= */
  t_focal = (epsilon - r.o.x)/r.d.x;
  focal_length = abs(t_inf - t_focal);

  /* ==== near_focus ========= */
  //Vector3D nearest_o = Vector3D(0,0,-5*focal_length);
  Vector3D nearest_o = Vector3D(0,0,elts.back().center - elts.back().radius - (1 + log(focal_length))*focal_length);
  cout << "close object distance:" << elts.back().center - elts.back().radius - (1 + log(focal_length))*focal_length << endl;
  r = Ray(nearest_o, (Vector3D(epsilon, 0, front_lens_z)-nearest_o).unit());
  trace.clear();
  trace.push_back(r.o);
  if (trace_backwards(r, &trace)) {
       if ((r.o.x !=0) && (r.d.x == 0)) {
          t_near = INF_D;
          near_focus = INF_D;
       } else {
          t_near = -r.o.x/r.d.x;
          near_focus = r.o.z + t_near*r.d.z;
       }
   } else {
       near_focus = -42;
   }
 
  cout << "[Lens] Infinity focus depth is " << infinity_focus << endl;
  cout << "[Lens] Close focus depth is " << near_focus << endl;
  cout << "[Lens] True focal length is " << focal_length << endl;

  cout << "Deliverable, sensor_depth, object_depth:" << endl;
  double step = (near_focus - infinity_focus)/100.0;
  cout << "sensor = [";
  for (double d = near_focus; d >= infinity_focus; d -= step) {
    cout << d << ",";
  }
  cout << "]" << endl << "obj = [";
  for (double d = near_focus; d >= infinity_focus; d -= step) {
    cout << -focus_depth(d) << ",";
  }
  cout << "]" << endl << "End Deliverable:" << endl;
}




bool Lens::trace(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the sensor out into the world.
  double prev_ior;
  prev_ior = 1.0;
  for (int i=0; i<elts.size(); i++) {
    LensElement elt = elts[i];
    if (!elt.pass_through(r, prev_ior)) return false;
    trace->push_back(r.o);        
   }
   return true;
}

bool Lens::trace_backwards(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the world backwards through 
  // the lens towards the sensor.
  double prev_ior;
  for (int i=elts.size()-1; i>=0; i--) {
    LensElement elt = elts[i];
    prev_ior = (i>0) ? elts[i-1].ior: 1.0;
    if (!elt.pass_through(r, prev_ior)) return false;
    trace->push_back(r.o);        
  }
  return true;
}

float Lens::focus_depth(float d) const {

  // Part 1 Task 2: Implement this. Should find the conjugate of a ray
  // starting from the sensor at depth d.
  double epsilon = ap_radius/150.0;
  LensElement back_lens = elts[0];
  double back_lens_z = back_lens.center-back_lens.radius;
  double t;
  Vector3D o = Vector3D(0,0,d);
  Ray r = Ray(Vector3D(0,0,d), Vector3D(epsilon, 0, back_lens_z-d).unit());
  vector<Vector3D> trace;
  trace.push_back(r.o);
  if (this->trace(r, &trace)) {
      if ((r.o.x !=0) && (r.d.x == 0)) {
         return INF_D;
      } else {
         t = -r.o.x/r.d.x;
         return r.o.z + t*r.d.z;
      }
  }
  return 0;
}

Vector3D Lens::back_lens_sample() const {

  // Part 1 Task 2: Implement this. Should return a point randomly sampled
  // on the back element of the lens (the element closest to the sensor)
  double radius = elts[0].aperture * 0.5;
  double z = elts[0].center-elts[0].radius;
  while (1) {
    double x = random_uniform()*2*radius - radius;
    double y = random_uniform()*2*radius - radius;
    if (x*x+y*y<=radius*radius)
      return Vector3D(x,y,z);
  }
}



/****** LensCamera functions ******/


LensCamera::LensCamera(): pt(NULL) {
  string path = string(__FILE__).substr(0,string(__FILE__).find_last_of('/')+1) + "../lenses/";
  static const vector<string> lens_files = {"dgauss.50mm.dat", "wide.22mm.dat", "telephoto.250mm.dat", "fisheye.10mm.dat"};
  for (string lens_file : lens_files)
    lenses.emplace_back(path + lens_file);

  mount_lens(0);
}


Ray LensCamera::generate_ray(double x, double y) const {

  Ray r = Ray(Vector3D(),Vector3D() );
  if (lens_ind >= 0) {

    // Part 1 Task 2: Implement this. It generates a ray from sensor pixel (x,y)
    // pointing toward the back element of the lens (use back_lens_sample) and traces
    // it through the Lens (using your "trace" function)
    double film_d = sqrt(24*24+36*36);
    double screen_d = sqrt(screenW*screenW + screenH*screenH);
    double film_w = film_d * screenW / screen_d;
    double film_h = film_d * screenH / screen_d;
    double sensor_depth = curr_lens().sensor_depth;
    Vector3D sensor_point(-(x-0.5)*film_w, -(y-0.5)*film_h, sensor_depth);
    Vector3D dir = (curr_lens().back_lens_sample() - sensor_point).unit();
    r = Ray(sensor_point,dir);
    vector<Vector3D> trace;
    trace.push_back(r.o);
    if (!curr_lens().trace(r, &trace)) {
      return Ray(sensor_point,Vector3D(0,0,1));
    }

    /***** end of your code ******/


    // This code converts the ray you traced through the lens into world coordinates.
    r.o = pos + c2w * r.o * scale;
    r.d = (c2w * r.d).unit();

  } else {

    // Generate ray for a pinhole camera. Same as in the previous assignment.
    x = 2*(x-.5); y = 2*(y-.5);
    r = Ray(pos,(c2w*Vector3D(x*tan(radians(hFov)*.5),y*tan(radians(vFov)*.5),-1)).unit());

  }

  r.min_t = nClip; r.max_t = fClip;
  return r;
}



void LensCamera::move_sensor(float delta) {
  if (lens_ind < 0) return;
  curr_lens().sensor_depth += delta;
  cout << "[LensCamera] Sensor plane moved to " << curr_lens().sensor_depth
       << ", focus now at " << lenses[lens_ind].focus_depth(lenses[lens_ind].sensor_depth) << endl;
  //double infinity_focus = curr_lens().infinity_focus;
  //double near_focus = curr_lens().near_focus;
  //double step = (near_focus - infinity_focus)/100.0;
  //for (double d = near_focus; d >= infinity_focus; d -= step) {
  //  cout << d << "," << lenses[lens_ind].focus_depth(d) << endl;
  //}
  //cout << "End Deliverable:" << endl;
}

void LensCamera::stop_down(float ratio) {
  float ap = curr_lens().ap_radius * ratio;
  if (ap > curr_lens().ap_original) ap = curr_lens().ap_original;
  curr_lens().ap_radius = ap;
  cout << "[LensCamera] Aperture is now " << curr_lens().ap_radius << "mm" << endl;
}

void LensCamera::mount_lens(int i) {
  lens_ind = i;
  if (i >= 0) {
    cout << "[LensCamera] Switched to lens #" << (i+1) 
         << " with focal length " << curr_lens().focal_length << "mm" << endl;
  } else {
    cout << "[LensCamera] Switched to pinhole camera" << endl;
  }
}



// A dummy function to demonstrate how to work with the image buffer.
// Calculates the average value of the green color channel in the image.
// You'll have to remember your 2D array indexing in order to take differences
// of neighboring pixels in a more sophisticated metric function.
static double mean_green(const ImageBuffer& ib) {
  double sum = 0;
  for (int i = 0; i < ib.w * ib.h; ++i) {
      sum += green_channel(ib.data[i]);
  }
  double mean = sum / (ib.w * ib.h);
  
  return mean;
}

double LensCamera::focus_metric(const ImageBuffer& ib) const {

  // Part 2 Task 1: Implement this. Design a metric to judge how "in-focus"
  // the image patch stored in the provided ImageBuffer is.
  Vector3D sum(0,0,0);
  Vector3D var(0,0,0);
  for (int i = 0; i < ib.w * ib.h; ++i) {
      sum.x += red_channel(ib.data[i]);
      sum.y += green_channel(ib.data[i]);
      sum.z += blue_channel(ib.data[i]);
  }
  Vector3D mean = sum / (ib.w * ib.h);
  for (int i = 0; i < ib.w * ib.h; ++i) {
      var.x += pow(red_channel(ib.data[i])-mean.x,2);
      var.y += pow(green_channel(ib.data[i])-mean.y,2);
      var.z += pow(blue_channel(ib.data[i])-mean.z,2);
  }
  var = var / (ib.w * ib.h);
  return (var.x+var.y+var.z)/3;
}


void LensCamera::autofocus() {


  // Part 2 Task 2: Implement this. Design a global search using your 
  // focus metric to set the sensor to be at the depth where the 
  // render cell is most "in focus". Provided code shows how to 
  // move the sensor, request a render of the cell, and evaluate the focus metric.

  // This call ensures that your pathtracer is rendering at high enough quality.
  // Increase samples per pixel to 16 and samples per light to 16.
  pt->bump_settings();

  // Example code. Nothing to do with your actual implementation except to 
  // demonstrate functionality.
  ImageBuffer ib;

  //cout << "[LensCamera] The mean green is " << focus_metric(ib) << endl;
  double sensor_pixel_size = sqrt(36*36 + 24*24) / sqrt(screenW*screenW + screenH*screenH);
  double zi_over_A = 2;
  double step = sensor_pixel_size * zi_over_A;
  cout << "[LensCamera] step:" << step << endl;
  double infinity_focus = curr_lens().infinity_focus;
  cout << "[LensCamera] inf:" << infinity_focus << endl;
  double near_focus = curr_lens().near_focus;
  cout << "[LensCamera] near:" << near_focus << endl;
  double d = infinity_focus;
  double max_var = -1;
  double best_d = 1;
  double metric;
  cout << "DeliveraBLE: d ; metric" << endl;
  cout << "d_m = [";
    while (d <= near_focus) {
    curr_lens().sensor_depth = d;
    pt->raytrace_cell(ib);
    metric = focus_metric(ib); 
    cout << "[" << d << "," << metric << "],";
    if (metric > max_var) {
        max_var = metric;
        best_d = d;
    }
    d += step;
  }
  cout << "]" << endl;
  curr_lens().sensor_depth = best_d;
  cout << "[LensCamera] best d:" << best_d << endl;
}





void LensCamera::dump_settings(string filename) {
  ofstream file(filename);
  file << hFov << " " << vFov << " " << ar << " " << nClip << " " << fClip << endl;
  for (int i = 0; i < 3; ++i)
    file << pos[i] << " ";
  for (int i = 0; i < 3; ++i)
    file << targetPos[i] << " ";
  file << endl;
  file << phi << " " << theta << " " << r << " " << minR << " " << maxR << endl;
  for (int i = 0; i < 9; ++i)
    file << c2w(i/3, i%3) << " ";
  file << endl;
  file << screenW << " " << screenH << " " << screenDist << endl;

  file << lens_ind << endl;
  for (Lens &lens : lenses) {
    file << lens.sensor_depth << " ";
  }
  file << endl;

  cout << "[LensCamera] Dumped settings to " << filename << endl;
}

void LensCamera::load_settings(string filename) {
  ifstream file(filename);

  file >> hFov >> vFov >> ar >> nClip >> fClip;
  for (int i = 0; i < 3; ++i)
    file >> pos[i];
  for (int i = 0; i < 3; ++i)
    file >> targetPos[i];
  file >> phi >> theta >> r >> minR >> maxR;
  for (int i = 0; i < 9; ++i)
    file >> c2w(i/3, i%3);
  file >> screenW >> screenH >> screenDist;

  file >> lens_ind;
  for (Lens &lens : lenses) {
    file >> lens.sensor_depth;
  }

  cout << "[LensCamera] Loaded settings from " << filename << endl;
}


} // namespace CGL

