#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, vector<double> maps_x, vector<double> maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2( (map_y-y),(map_x-x) );

	double angle = abs(theta-heading);

	if(angle > pi()/4)
	{
		closestWaypoint++;
	}

	return closestWaypoint;

}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, vector<double> maps_x, vector<double> maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, vector<double> maps_s, vector<double> maps_x, vector<double> maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

// Get lane Id given Frenet d coordinates
int getLane (double d){
    
    int laneId = -1;
    
    if( d > 0 && d < 4)
        laneId = 0;
    else if (d > 4 && d < 8)
        laneId = 1;
    else if (d > 8 && d < 12)
        laneId = 2;
    
    return laneId;
}

// Check for cars ahead and behind given a lane
vector<double> contiguousLaneBehaviour(vector<vector<double >> sensor_fusion, int lane, int prev_size, double car_s, double car_d, double end_path_s, double end_path_d, int forward_offset){
    
    double status_car_ahead = 0; // 0: no car ahead, 1: car ahead in range
    double status_car_behind = 0; // 0: no car behind, 1: car behind in range
    double dist_car_ahead = 10000;
    double dist_car_behind = 10000;
    double v_car_ahead = 0;
    double v_car_behind = 0;
    double id_car_ahead = -1;
    double id_car_behind = -1;

    double car_end_s = prev_size > 0 ? end_path_s : car_s;
    
    for(int i = 0; i < sensor_fusion.size(); i++){
    
        // [id, x, y, v_x, v_y, s, d]
        vector<double> i_car = sensor_fusion[i];
        double i_car_d = i_car[6];
        
        if (i_car_d > (2 + 4 * lane - 2) && i_car_d < (2 + 4 * lane + 2)){
            
            // ith car is in the same lane as in the parameter lane
            double i_car_v_x = i_car[3];
            double i_car_v_y = i_car[4];
            double i_car_v = sqrt(i_car_v_x * i_car_v_x + i_car_v_y * i_car_v_y); // ith car speed
            double i_car_s = i_car[5];
            double i_car_end_s = i_car_s + prev_size * 0.02 * i_car_v; // predicted ith car end s location
            
            if (i_car_end_s > car_end_s && i_car_end_s < (car_end_s + 30 + forward_offset)){ // Add offset so as to avoid lane changes that are not really useful for overtaking, i.e. when the car ahead in the new lane is slightly ahead than the car ahead in the current lane
                // The i_th car is ahead and in the path planning horizon
                status_car_ahead = 1;
                double i_dist_car_ahead = i_car_end_s - car_end_s;
                if (i_dist_car_ahead < dist_car_ahead){
                    dist_car_ahead = i_dist_car_ahead;
                    v_car_ahead = i_car_v;
                    id_car_ahead = i_car[0];
                }
            }
            else if (i_car_end_s < car_end_s && car_end_s < (i_car_end_s + 30)){
                // The i_th car is behind and in the path planning horizon
                status_car_behind = 1;
                double i_dist_car_behind = car_end_s - i_car_end_s;
                if (i_dist_car_behind < dist_car_behind){
                    dist_car_behind = i_dist_car_behind;
                    v_car_behind = i_car_v;
                    id_car_behind = i_car[0];
                }
            }
        }
    }
    
    vector<double> lane_behaviour = {status_car_ahead, status_car_behind, dist_car_ahead, dist_car_behind, v_car_ahead, v_car_behind, id_car_ahead, id_car_behind};
    return lane_behaviour;
    
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }

  // Lane the car is currently on
  int lane = 1;
    
  // Target speed
  double ref_speed = 0; // mph
  double max_speed = 49; // mph
    
  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,&lane,&ref_speed,&max_speed](uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);
        
        string event = j[0].get<string>();
        
        if (event == "telemetry") {
          // j[1] is the data JSON object
          
        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];

          	// Previous path data given to the Planner
          	auto previous_path_x = j[1]["previous_path_x"];
          	auto previous_path_y = j[1]["previous_path_y"];
          	// Previous path's end s and d values 
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

            int prev_size = previous_path_x.size();
            
          	// Define a path made up of (x,y) points that the car will visit sequentially every .02 seconds
            
            // Update car reference. Too smooth trajectories, use end_path of previous path as the starting reference point in the current iteration
            if(prev_size > 0){
                car_s = end_path_s;
            }
            
            // Update target lane and target_speed with sensor fusion data
            bool too_close = false;
            
            // Iterate through all of sensor fusion cars
            for(int i = 0; i < sensor_fusion.size(); i++){
                // ith car is in current target lane
                double i_d = sensor_fusion[i][6];
                if(i_d < (2 + 4 * lane + 2) && i_d > (2 + 4 * lane - 2)){
                    double i_vx = sensor_fusion[i][3];
                    double i_vy = sensor_fusion[i][4];
                    double i_v = sqrt(i_vx*i_vx + i_vy * i_vy);
                    double i_s = sensor_fusion[i][5];
                    
                    // Predict ith car s value, using the time interval (0.02 sec) times the the waypoints in the previous path
                    double pred_i_s = i_s + (double) prev_size * 0.02 * i_v;
                    
                    // If predicted ith car s value is too close to the end waypoint of the previous path, then slow down or change lanes
                    
                    if((pred_i_s > car_s) && (pred_i_s - car_s < 30)){
                        // TODO: use FSM cost function logic to decide to wether prepare lane change left or right, or just slow down
                        too_close = true;
                    }
                }
            }
            
            if(too_close){
                ref_speed -= 0.224;
            }
            else if (ref_speed < max_speed){
                ref_speed += 0.224;
            }
            
            // Find
          	
            // As in walkthrough video, create a auxiliary list of sparse (x,y) waypoints (aux_wps), evenly spaced at 30m
            // We will later use this auxiliary waypoints to interporlate and create more waypoints so as to smooth the trajectory
            vector<double> pts_x;
            vector<double> pts_y;
            
            double ref_x = car_x;
            double ref_y = car_y;
            double ref_yaw = deg2rad(car_yaw);
            
            // If previous path is almost empty, use car coordinates and heading as a starting reference to fill auxiliary list
            if(prev_size < 2){
                
                // Use points tangent to car heading and add them to the auxiliary list
                // TODO: ¿Why not multiply cos/sin by car speed and time?
                double prev_car_x = car_x - cos(car_yaw);
                double prev_car_y = car_y - sin(car_yaw);
                
                pts_x.push_back(prev_car_x);
                pts_x.push_back(car_x);
                
                pts_y.push_back(prev_car_y);
                pts_y.push_back(car_y);
                
            }
            else{
            
                // Redefine reference state using the latest 2 points in the previous path
                ref_x = previous_path_x[prev_size-1];
                ref_y = previous_path_y[prev_size-1];
                
                double ref_x_prev = previous_path_x[prev_size-2];
                double ref_y_prev = previous_path_y[prev_size-2];
                ref_yaw = atan2(ref_y-ref_y_prev, ref_x-ref_x_prev);
                
                pts_x.push_back(ref_x_prev);
                pts_x.push_back(ref_x);
                
                pts_y.push_back(ref_y_prev);
                pts_y.push_back(ref_y);
            }
            
            // Using Frenet coordinates, add evenly spaced points (30m) ahead of the starting reference
            vector<double> next_wp_0 = getXY(car_s + 30, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp_1 = getXY(car_s + 60, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);
            vector<double> next_wp_2 = getXY(car_s + 90, (2 + 4 * lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

            pts_x.push_back(next_wp_0[0]);
            pts_x.push_back(next_wp_1[0]);
            pts_x.push_back(next_wp_2[0]);
            
            pts_y.push_back(next_wp_0[1]);
            pts_y.push_back(next_wp_1[1]);
            pts_y.push_back(next_wp_2[1]);
            
            json msgJson;
            
            // Define the points that will be used by the trajectory planner
            vector<double> next_x_vals;
            vector<double> next_y_vals;
            
            // Interpolate the auxiliary waypoints (aux_wps) using a spline, but before shift the aux_wps to vehicle coordinates (x-axis is aligned with the car heading)
            for(int i = 0; i < pts_x.size(); i++){
                double shift_x = pts_x[i]-ref_x;
                double shift_y = pts_y[i]-ref_y;
                
                cout << "pts_x prev " << i << " " << pts_x[i] << endl;
                cout << "pts_y prev " << i << " " << pts_y[i] << endl;
                
                pts_x[i] = shift_x * cos(-ref_yaw) - shift_y * sin(-ref_yaw);
                pts_y[i] = shift_x * sin(-ref_yaw) + shift_y * cos(-ref_yaw);
                
                cout << "pts_x " << i << " " << pts_x[i] << endl;
                cout << "pts_y " << i << " " << pts_y[i] << endl;
            }
            
            // Create spline and set (x,y) points to the spline
            tk::spline sp;
            
            sp.set_points(pts_x,pts_y);
            
            // Fill first the previous path points from the last iterarion, i.e., the remaining points not used by the simulator
            for(int i = 0; i < previous_path_x.size(); i++){
                next_x_vals.push_back(previous_path_x[i]);
                next_y_vals.push_back(previous_path_y[i]);
            }
            
            // Use the method in the walkthrough in which the spline points projected along the x-axis are splitted into equal sizes so as to travel at the target speed
            double sp_target_x = 30.0;
            double sp_target_y = sp(sp_target_x);
            double sp_target_dist = sqrt((sp_target_x)*(sp_target_x) + (sp_target_y)*(sp_target_y));
            
            double x_add_on = 0;
            
            // Fill up the remaining waypoints with of the trajectory planner
            for(int i = 1; i <= 50 - previous_path_x.size(); i++){
            
                // N is the number of splits, which is calculated using the linearized distance (sp_target_dist), the time interval in which a waypoint is visited (0.02) and the target speed (using the 2.24 to convert from mph to mps).
                // TODO: update target speed in this loop as well
                double N = (sp_target_dist/(0.02 * ref_speed / 2.24));
                double x_point = x_add_on + sp_target_x / N;
                double y_point = sp(x_point);
                
                x_add_on = x_point;
                
                // Rotate back to normal coordinates
                double x_ref = x_point;
                double y_ref = y_point;
                
                x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
                y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));
                
                x_point += ref_x;
                y_point += ref_y;
                
                next_x_vals.push_back(x_point);
                next_y_vals.push_back(y_point);
            }
            
            // END
            
            msgJson["next_x"] = next_x_vals;
          	msgJson["next_y"] = next_y_vals;

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
          
        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
















































































