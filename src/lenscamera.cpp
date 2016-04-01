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
  

  return true;
}
bool LensElement::intersect(const Ray &r, Vector3D *hit_p) const {
  // Part 1 Task 1: Implement this. It intersects r with this spherical lens elemnent 
  // (or aperture diaphragm). You'll want to reuse some sphere intersect code.


  return true;
  
}
bool LensElement::refract(Ray& r, const Vector3D& hit_p, const double& prev_ior) const {
  // Part 1 Task 1: Implement this. It refracts the Ray r with this lens element or 
  // does nothing at the aperture element.
  // You'll want to consult your refract function from the previous assignment.



  return true;

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



  cout << "[Lens] Infinity focus depth is " << infinity_focus << endl;
  cout << "[Lens] Close focus depth is " << near_focus << endl;
  cout << "[Lens] True focal length is " << focal_length << endl;
}




bool Lens::trace(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the sensor out into the world.



  return true;
}

bool Lens::trace_backwards(Ray &r, std::vector<Vector3D> *trace) const {
  // Part 1 Task 1: Implement this. It traces a ray from the world backwards through 
  // the lens towards the sensor.



  return true;
}

float Lens::focus_depth(float d) const {

  // Part 1 Task 2: Implement this. Should find the conjugate of a ray
  // starting from the sensor at depth d.

  return 0;
}

Vector3D Lens::back_lens_sample() const {

  // Part 1 Task 2: Implement this. Should return a point randomly sampled
  // on the back element of the lens (the element closest to the sensor)

  return Vector3D();

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

  return mean_green(ib); //  A meaningless standin
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
  curr_lens().sensor_depth += 1;
  pt->raytrace_cell(ib);
  cout << "[LensCamera] The mean green is " << focus_metric(ib) << endl;


  
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

