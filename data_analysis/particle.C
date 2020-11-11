/*
Class to store particle 4-momenta and calculate useful quantities
*/
#include "math.h"

class particle{

    private:
        float m_x_momenta;
        float m_y_momenta;
        float m_z_momenta;
        float m_trans_momenta;
        float m_mass;
        float m_energy;
        float m_pseudorapidity;
        float m_rapidity;
        float m_theta;
        float m_phi;
        int m_pdgid;
        int m_set_method;
        std::string m_type;
        std::vector<particle> m_particles;

    public:
        particle() : m_x_mometa{0}, m_y_mometa{0}, m_z_mometa{0}, m_trans_momenta{0}, m_mass{0}, m_energy{0}, m_pseudorapidity{0}, m_rapidity{0}, m_theta{0}, m_phi{0}, m_pdgid{0}, m_set_method{0}, m_type{"None"}, m_particles{0} {}
        ~particle()

        void set_PtEtaPhiM(const float&, const float&, const float&, const float&);
        void set_EPxPyPz(const float&, const float&, const float&, const float&);
        void set_pdgid(const int &);
        void set_type(const std::string&);
        float x_momenta();
        float y_momenta();
        float z_momenta();
        float total_momenta();
        float energy();
        float mass();
        float pseudorapidity();
        float rapidty();
        float theta();
        float phi();
        int pdgid();
        std::string type();
        std::vector<particle> particles();
        float deltaR(const particle&);
        float centrality(const particle&, const particle&);
        void append(const particle&);
        void sort_by_pT();

        particle& operator[](const int&);
        particle& operator+(const particle&);
        particle& operator-(const particle&);
        particle& operator=(const particle&);

}

void particle::set_PtEtaPhiM(const float& trans_momenta, const float& pseudorapidity, const float& phi, const float& mass){
    m_trans_momenta = trans_momenta; m_pseudorapidity = pseudoRapidity; m_phi = phi, m_mass = mass; m_set_method = 1;
    m_particles.push_back(&this)
}
void particle::set_EPxPyPZ(const float& energy, const float& x_momenta, const float& y_momenta, const float& z_momenta){
    m_energy = energy; m_x_momenta = x_momenta; m_y_momenta = y_momenta; m_z_momenta = z_momenta; m_set_method = 2;
    m_particles.push_back(&this)
}
void particle::set_pdgid(const int& pdgid){m_pdgid = pdgid; }
float particle::x_momenta(){
    if(m_set_method == 1){ return m_trans_momenta *  cos(m_phi); }
    if(m_set_method == 2){ return m_x_momenta; }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::y_momenta(){
    if(m_set_method == 1){ return m_trans_momenta *  sin(m_phi); }
    if(m_set_method == 2){ return m_y_momenta; }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::z_momenta(){
    if(m_set_method == 1){ return m_trans_momenta *  sinh(m_pseudorapidity); }
    if(m_set_method == 2){ return m_z_momenta; }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::total_momenta(){
    if(m_set_method == 1){ return m_trans_momenta *  cosh(m_pseudorapidity); }
    if(m_set_method == 2){ return  sqrt(pow(m_x_momenta, 2) + pow(m_y_momenta, 2), pow(m_z_momenta, 2)); }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::energy(){
    if(m_set_method == 1){ return  sqrt(pow(m_mass, 2) + pow(m_trans_momenta *  cosh(m_phi), 2)); }
    if(m_set_method == 2){ return m_energy; }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::trans_momenta(){
    if(m_set_method == 1){ return m_trans_momenta; }
    if(m_set_method == 2){ return  sqrt(pow(m_x_momenta, 2) + pow(m_y_momenta, 2)); }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::pseudorapidity(){   
    if(m_set_method == 1){ return m_pseudorapidity; }
    if(m_set_method == 2){ return  atanh(m_z_momenta /  sqrt(pow(m_x_momenta, 2) + pow(m_y_momenta, 2), pow(m_z_momenta, 2))); }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::rapidty(){
    if(m_set_method == 1){ 
        numerator =  sqrt(pow(m_mass, 2) + pow(m_trans_momenta *  cosh(m_pseudorapidity), 2)) + m_trans_momenta * m_pseudorapidity;
        denominator =  sqrt(pow(m_mass, 2) + pow(m_pseudorapidity, 2));
        return  log(numerator / denominator);
    }
    if(m_set_method == 2){ return  log((m_energy + m_z_momenta) / (m_energy - m_z_momenta)); }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::phi(){
    if(m_set_method == 1){ return m_phi; }
    if(m_set_method == 2){ return  acos(m_x_momenta /  sqrt(pow(m_x_momenta, 2) + pow(m_y_label, 2))); }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
float particle::theta(){
    if(m_set_method == 1){ return 2 *  atan( exp(-m_pseudorapidity)); }
    if(m_set_method == 2){ return 2 *  atan( exp(- atanh(m_z_momenta /  sqrt(pow(m_x_momenta, 2) + pow(m_y_momenta, 2), pow(m_z_momenta, 2))))); }
    else{ std::cout << "ERROR: No set method was called for this particle" << std::endl; return 0; }
}
int particle::pdgid(){
    return m_pdgid; 
}
std::vector<particle> particle::particles(){
    return m_particles; 
}
std::string particle::type(){
    return m_type;
}
float particle::deltaR(const particle& particle2){
    return  sqrt(pow(this->pseudorapidity() - particle2.pseudorapidity(), 2) + pow(this->phi() - particle2.phi(), 2)); 
}
float particle::centrality(const particle& particle2, const particle& particle3){
    return (this->pseudorapidity() - ((particle2.pseudorapidity() + particle3.pseudorapidity()) / 2)) / (particle2.pseudorapidity() - particle3.pseudorapidity()) 
}
void particle::append(const particle& particle2){ 
    m_particles.push_back(particle2); 
}
void particle::sort_by_pT(){
    bool pTComp(particle lhs, particle rhs){ 
    return lhs.trans_momenta() > rhs.trans_momenta();
    }
    std::sort(m_particles.begin(), particles.end(), pTcomp);
    if(m_particles[0].m_set_method = 1){this->set_PtEtaPhiM(m_particles[0].m_trans_momenta, m_particles[0].m_pseudorapidity, m_particles[0].m_phi, m_particles[0].m_mass); }
    if(m_particles[0].m_set_method = 2){this->set_EPxPyPz(m_particles[0].m_energy, m_particles[0].m_x_momenta, m_particles[0].m_y_momenta, m_particles[0].m_z_momenta); }
}



particle& particle::operator[](int index){
    if(abs(index) < m_particles.size()){
        if(index > 0){ return m_particles[index]; }
        if(index < 0){ return m_particles[m_particles.size()-abs(index)-1]; }
    }
    else{
        std::cout << "Error: Index " << index << " is out of bounds for a vector or particles of size " << m_particles.size() << std::endl;
        return particle();
    }
}

particle& particle::operator+(const particle& particle2){
    particle new_particle; 
    new_particle.m_x_momenta = this->x_momenta() + particle2.x_momenta();
    new_particle.m_y_momenta = this->y_momenta() + particle2.y_momenta();
    new_particle.m_z_momenta = this->z_momenta() + particle2.z_momenta();
    new_particle.energy = this->energy() + particle2.energy();
    return new_particle;
}

particle& particle::operator-(const particle& particle2){
    particle new_particle; 
    new_particle.m_x_momenta = this->x_momenta() - particle2.x_momenta();
    new_particle.m_y_momenta = this->y_momenta() - particle2.y_momenta();
    new_particle.m_z_momenta = this->z_momenta() - particle2.z_momenta();
    new_particle.energy = this->energy() - particle2.energy();
    return new_particle;
}

particle& particle::operator=(const particle& particle2){
    if(this != &particle2){
        m_x_momenta = particle2.x_momenta();
        m_y_momenta = particle2.y_momenta();
        m_z_momenta = particle2.z_momenta();
        m_energy = particle2.energy();
        m_particles.clear();
        m_particles = particle2.particles();
        }
    return *this
}
