/*
This program simulates the three-body problem using the Leapfrog integration method.
It visualizes the motion of three celestial bodies under mutual gravitational attraction
in a 2D space. The program includes:

1. Leapfrog integration for numerical stability and efficiency.
2. Runge-Kutta method for acceleration calculation (not currently used).
   - I initially tried using the Runge-Kutta method because I had learned about it before.
     However, since it isn't symplectic, it's not well-suited for chaotic simulations like the three-body problem.
     This was an interesting insight I discovered while working on the project.
3. A "Body" class to represent each celestial body, including position, velocity, mass, and rendering.
4. Trails for each body to visualize their paths over time.
5. A scaling factor to adjust the simulation's size and gravitational constant.

The simulation runs in a graphical window using the SFML library, where the bodies
move according to the laws of gravity, and their trails are drawn to show their trajectories.
*/

#include <iostream>
#include <SFML/Graphics.hpp>
#include <vector>
#include <cmath>
#include <deque>

// Scaling factor for the simulation
double scalingfactor = 3;
// Gravitational constant, scaled for the simulation
double G = 1 * pow(scalingfactor, 2);

// Function to calculate acceleration using the Runge-Kutta method
sf::Vector2<double> RungeKutta_K_N(sf::Vector2<double> pos1, sf::Vector2<double> pos2, sf::Vector2<double> pos3, double m1, double m2) {
    // Calculate distances between bodies
    float distance1 = sqrt(pow((pos1.x - pos2.x), 2) + pow((pos1.y - pos2.y), 2));
    float distance2 = sqrt(pow((pos1.x - pos3.x), 2) + pow((pos1.y - pos3.y), 2));

    // Calculate acceleration due to gravitational forces
    sf::Vector2<double> kn_vn = -G * ((m1) * (pos1 - pos2) / pow(distance1, 3) + (m2) * (pos1 - pos3) / pow(distance2, 3));
    return kn_vn;
}

// Function to calculate new positions using the Runge-Kutta method
std::vector<sf::Vector2<double>> calculateNewPosition(double dt, sf::Vector2<double> vel1, sf::Vector2<double> vel2, sf::Vector2<double> vel3, sf::Vector2<double> pos1, sf::Vector2<double> pos2, sf::Vector2<double> pos3) {
    sf::Vector2<double> k1_r1 = vel1;
    sf::Vector2<double> k1_r2 = vel2;
    sf::Vector2<double> k1_r3 = vel3;

    sf::Vector2<double> k2_r1 = (k1_r1 * (dt / 2) + vel1);
    sf::Vector2<double> k2_r2 = (k1_r2 * (dt / 2) + vel2);
    sf::Vector2<double> k2_r3 = (k1_r3 * (dt / 2) + vel3);

    sf::Vector2<double> k3_r1 = (k2_r1 * (dt / 2) + vel1);
    sf::Vector2<double> k3_r2 = (k2_r2 * (dt / 2) + vel2);
    sf::Vector2<double> k3_r3 = (k2_r3 * (dt / 2) + vel3);

    sf::Vector2<double> k4_r1 = (k3_r1 * dt + vel1);
    sf::Vector2<double> k4_r2 = (k3_r2 * dt + vel2);
    sf::Vector2<double> k4_r3 = (k3_r3 * dt + vel3);

    return {
        pos1 + (dt / 6) * (k1_r1 + 2.0 * k2_r1 + 2.0 * k3_r1 + k4_r1),
        pos2 + (dt / 6) * (k1_r2 + 2.0 * k2_r2 + 2.0 * k3_r2 + k4_r2),
        pos3 + (dt / 6) * (k1_r3 + 2.0 * k2_r3 + 2.0 * k3_r3 + k4_r3)
    };
}

// Function to calculate new velocities using the Runge-Kutta method
std::vector<sf::Vector2<double>> calculateNewVelocity(double dt, sf::Vector2<double> pos1, sf::Vector2<double> pos2, sf::Vector2<double> pos3, double m1, double m2, double m3, sf::Vector2<double> vel1, sf::Vector2<double> vel2, sf::Vector2<double> vel3) {
    sf::Vector2<double> k1_v1 = RungeKutta_K_N(pos1, pos2, pos3, m2, m3);
    sf::Vector2<double> k1_v2 = RungeKutta_K_N(pos2, pos1, pos3, m1, m3);
    sf::Vector2<double> k1_v3 = RungeKutta_K_N(pos3, pos1, pos2, m1, m2);

    sf::Vector2<double> k2_v1 = RungeKutta_K_N(pos1 + (dt / 2) * k1_v1, pos2 + (dt / 2) * k1_v2, pos3 + (dt / 2) * k1_v3, m2, m3);
    sf::Vector2<double> k2_v2 = RungeKutta_K_N(pos2 + (dt / 2) * k1_v2, pos1 + (dt / 2) * k1_v1, pos3 + (dt / 2) * k1_v3, m1, m3);
    sf::Vector2<double> k2_v3 = RungeKutta_K_N(pos3 + (dt / 2) * k1_v3, pos1 + (dt / 2) * k1_v1, pos2 + (dt / 2) * k1_v2, m1, m2);

    sf::Vector2<double> k3_v1 = RungeKutta_K_N(pos1 + (dt / 2) * k2_v1, pos2 + (dt / 2) * k2_v2, pos3 + (dt / 2) * k2_v3, m2, m3);
    sf::Vector2<double> k3_v2 = RungeKutta_K_N(pos2 + (dt / 2) * k2_v2, pos1 + (dt / 2) * k2_v1, pos3 + (dt / 2) * k2_v3, m1, m3);
    sf::Vector2<double> k3_v3 = RungeKutta_K_N(pos3 + (dt / 2) * k2_v3, pos1 + (dt / 2) * k2_v1, pos2 + (dt / 2) * k2_v2, m1, m2);

    sf::Vector2<double> k4_v1 = RungeKutta_K_N(pos1 + dt * k3_v1, pos2 + dt * k3_v2, pos3 + dt * k3_v3, m2, m3);
    sf::Vector2<double> k4_v2 = RungeKutta_K_N(pos2 + dt * k3_v2, pos1 + dt * k3_v1, pos3 + dt * k3_v3, m1, m3);
    sf::Vector2<double> k4_v3 = RungeKutta_K_N(pos3 + dt * k3_v3, pos1 + dt * k3_v1, pos2 + dt * k3_v2, m1, m2);

    return {
        vel1 + (dt / 6) * (k1_v1 + 2.0 * k2_v1 + 2.0 * k3_v1 + k4_v1),
        vel2 + (dt / 6) * (k1_v2 + 2.0 * k2_v2 + 2.0 * k3_v2 + k4_v2),
        vel3 + (dt / 6) * (k1_v3 + 2.0 * k2_v3 + 2.0 * k3_v3 + k4_v3)
    };
}

// Function to calculate acceleration for the Leapfrog method
sf::Vector2<double> Leapfrog_Acceleration(sf::Vector2<double> pos1, sf::Vector2<double> pos2, sf::Vector2<double> pos3, double m1, double m2) {
    // Calculate distances between bodies
    float distance1 = sqrt(pow((pos1.x - pos2.x), 2) + pow((pos1.y - pos2.y), 2));
    float distance2 = sqrt(pow((pos1.x - pos3.x), 2) + pow((pos1.y - pos3.y), 2));

    // Calculate acceleration due to gravitational forces
    sf::Vector2<double> a = -G * ((m1) * (pos1 - pos2) / pow(distance1, 3) + (m2) * (pos1 - pos3) / pow(distance2, 3));
    return a;
}

// Function to calculate half-step velocity for the Leapfrog method
sf::Vector2<double> Leapfrog_VelocityHalf(double dt, sf::Vector2<double> a, sf::Vector2<double> vel) {
    return vel + dt / 2 * a;
}

// Function to calculate new position for the Leapfrog method
sf::Vector2<double> Leapfrog_Position(double dt, sf::Vector2<double> velhalf, sf::Vector2<double> pos) {
    return pos + dt * velhalf;
}

// Function to calculate new positions and half-step velocities using the Leapfrog method
std::vector<sf::Vector2<double>> LeapfrogNewPos(double dt, sf::Vector2<double> vel1, sf::Vector2<double> vel2, sf::Vector2<double> vel3, sf::Vector2<double> pos1, sf::Vector2<double> pos2, sf::Vector2<double> pos3, double m1, double m2, double m3) {
    // Calculate accelerations
    sf::Vector2<double> a1 = Leapfrog_Acceleration(pos1, pos2, pos3, m2, m3);
    sf::Vector2<double> a2 = Leapfrog_Acceleration(pos2, pos1, pos3, m1, m3);
    sf::Vector2<double> a3 = Leapfrog_Acceleration(pos3, pos1, pos2, m1, m2);

    // Calculate half-step velocities
    sf::Vector2<double> v1_05 = Leapfrog_VelocityHalf(dt, a1, vel1);
    sf::Vector2<double> v2_05 = Leapfrog_VelocityHalf(dt, a2, vel2);
    sf::Vector2<double> v3_05 = Leapfrog_VelocityHalf(dt, a3, vel3);

    // Calculate new positions
    sf::Vector2<double> newpos1 = Leapfrog_Position(dt, v1_05, pos1);
    sf::Vector2<double> newpos2 = Leapfrog_Position(dt, v2_05, pos2);
    sf::Vector2<double> newpos3 = Leapfrog_Position(dt, v3_05, pos3);

    return {newpos1, newpos2, newpos3, v1_05, v2_05, v3_05};
}

// Function to calculate new velocities using the Leapfrog method
std::vector<sf::Vector2<double>> LeapfrogNewVel(double dt, sf::Vector2<double> pos1, sf::Vector2<double> pos2, sf::Vector2<double> pos3, double m1, double m2, double m3, sf::Vector2<double> vel1, sf::Vector2<double> vel2, sf::Vector2<double> vel3) {
    // Calculate new accelerations
    sf::Vector2<double> newa1 = Leapfrog_Acceleration(pos1, pos2, pos3, m2, m3);
    sf::Vector2<double> newa2 = Leapfrog_Acceleration(pos2, pos1, pos3, m1, m3);
    sf::Vector2<double> newa3 = Leapfrog_Acceleration(pos3, pos1, pos2, m1, m2);

    // Calculate full-step velocities
    sf::Vector2<double> v1 = Leapfrog_VelocityHalf(dt, newa1, vel1);
    sf::Vector2<double> v2 = Leapfrog_VelocityHalf(dt, newa2, vel2);
    sf::Vector2<double> v3 = Leapfrog_VelocityHalf(dt, newa3, vel3);
    return {v1, v2, v3};
}

// Class representing a celestial body
class Body {
private:
    sf::Vector2<double> velocity; // Velocity of the body
    sf::Vector2<double> position; // Position of the body
    double mass;                  // Mass of the body
    sf::CircleShape shape;        // Shape for rendering the body

public:
    // Constructor to initialize the body
    Body(double velX, double velY, double posX, double posY, double m, double r, sf::Color Color) {
        velocity = sf::Vector2<double>{velX, velY};
        position = sf::Vector2<double>{posX, posY};
        mass = m;
        shape.setRadius(r);
        shape.setPosition(400 + posX * 100, 400 + posY * 100);
        shape.setFillColor(Color);
    }

    // Set the position of the body
    void setPosition(sf::Vector2<double> pos) {
        position = pos;
        shape.setPosition((400 + position.x * 100), (400 + position.y * 100));
    }

    // Set the velocity of the body
    void setVelocity(sf::Vector2<double> vel1) {
        velocity = vel1;
    }

    // Get the position of the body
    sf::Vector2<double> getPosition() {
        return position;
    }

    // Get the position of the body in window coordinates
    sf::Vector2<double> getWindowPosition() {
        return {400 + position.x * 100, 400 + position.y * 100};
    }

    // Get the velocity of the body
    sf::Vector2<double> getVelocity() {
        return velocity;
    }

    // Get the mass of the body
    double getMass() {
        return mass;
    }

    // Draw the body on the window
    void draw(sf::RenderWindow &window) {
        window.draw(shape);
    }

    // Get the radius of the body
    double getRadius() {
        return shape.getRadius();
    }
};

int main() {
    sf::RenderWindow window(sf::VideoMode({800, 800}), "3-Body-Sim");

    // Time step for the simulation
    double h = 0.001f;

    // Initial positions and velocities for the figure-8 solution
    sf::Vector2<double> r1 = {scalingfactor * -0.97000436, scalingfactor * 0.24308753};
    sf::Vector2<double> r2 = {scalingfactor * 0.97000436, scalingfactor * -0.24308753};
    sf::Vector2<double> r3 = {scalingfactor * 0, scalingfactor * 0};

    sf::Vector2<double> v1 = {sqrt(scalingfactor) * 0.466203685, sqrt(scalingfactor) * 0.43236573};
    sf::Vector2<double> v2 = {sqrt(scalingfactor) * 0.466203685, sqrt(scalingfactor) * 0.43236573};
    sf::Vector2<double> v3 = {sqrt(scalingfactor) * -0.93240737, sqrt(scalingfactor) * -0.86473146};

    // Create the celestial bodies
    Body body1(v1.x, v1.y, r1.x, r1.y, 1, 20., sf::Color::Green);
    Body body2(v2.x, v2.y, r2.x, r2.y, 1, 20., sf::Color::Red);
    Body body3(v3.x, v3.y, r3.x, r3.y, 1, 20., sf::Color::Blue);

    // Deques to store the trails of the bodies
    std::deque<sf::Vertex> trail1;
    std::deque<sf::Vertex> trail2;
    std::deque<sf::Vertex> trail3;

    while (window.isOpen()) {
        // Calculate new positions and velocities using the Runge-Kutta method
        /*
        std::vector<sf::Vector2<double>> v = calculateNewVelocity(h, body1.getPosition(), body2.getPosition(), body3.getPosition(), body1.getMass(), body2.getMass(), body3.getMass(), body1.getVelocity(), body2.getVelocity(), body3.getVelocity());
        std::vector<sf::Vector2<double>> r = calculateNewPosition(h, body1.getVelocity(), body2.getVelocity(), body3.getVelocity(), body1.getPosition(), body2.getPosition(), body3.getPosition());
        */

        // Calculate new positions and velocities using the Leapfrog method
        
        std::vector<sf::Vector2<double>> r = LeapfrogNewPos(h, body1.getVelocity(), body2.getVelocity(), body3.getVelocity(), body1.getPosition(), body2.getPosition(), body3.getPosition(), body1.getMass(), body2.getMass(), body3.getMass());
        std::vector<sf::Vector2<double>> v = LeapfrogNewVel(h, r[0], r[1], r[2], body1.getMass(), body2.getMass(), body3.getMass(), r[3], r[4], r[5]);
        
        // Update the positions and velocities of the bodies
        body1.setPosition(r[0]);
        body1.setVelocity(v[0]);

        body2.setPosition(r[1]);
        body2.setVelocity(v[1]);

        body3.setPosition(r[2]);
        body3.setVelocity(v[2]);

        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            } else if (event.type == sf::Event::Resized) {
                unsigned int newWidth = event.size.width;
                unsigned int newHeight = event.size.height;

                sf::View view = window.getView();
                view.setSize(static_cast<float>(newWidth), static_cast<float>(newHeight));
                window.setView(view);
            }
        }

        window.clear();

        // Update trails with the new positions
        trail1.push_back(sf::Vertex(sf::Vector2f(body1.getWindowPosition().x + body1.getRadius(), body1.getWindowPosition().y + body1.getRadius()), sf::Color::Green));
        trail2.push_back(sf::Vertex(sf::Vector2f(body2.getWindowPosition().x + body2.getRadius(), body2.getWindowPosition().y + body2.getRadius()), sf::Color::Red));
        trail3.push_back(sf::Vertex(sf::Vector2f(body3.getWindowPosition().x + body3.getRadius(), body3.getWindowPosition().y + body3.getRadius()), sf::Color::Blue));

        // Limit trail size to avoid performance issues
        if (trail1.size() > 1000) trail1.pop_front();
        if (trail2.size() > 1000) trail2.pop_front();
        if (trail3.size() > 1000) trail3.pop_front();

        // Convert trails to vertex arrays for rendering
        sf::VertexArray trail1Array(sf::LineStrip, trail1.size());
        sf::VertexArray trail2Array(sf::LineStrip, trail2.size());
        sf::VertexArray trail3Array(sf::LineStrip, trail3.size());

        for (size_t i = 0; i < trail1.size(); ++i) trail1Array[i] = trail1[i];
        for (size_t i = 0; i < trail2.size(); ++i) trail2Array[i] = trail2[i];
        for (size_t i = 0; i < trail3.size(); ++i) trail3Array[i] = trail3[i];

        // Draw trails and bodies
        window.draw(trail1Array);
        window.draw(trail2Array);
        window.draw(trail3Array);

        body1.draw(window);
        body2.draw(window);
        body3.draw(window);

        window.display();
    }

    return 0;
}
