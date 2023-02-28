%{
Indirect trust (Dempster-Shaffer (D-S) theory) calculation 
using Evaluation Error, Similarity Parameters, Recommendation Credibiltiy, Belief Function
and Mass Functions

The following code is just a Snippet without any inputs, It has to be integrated into the
larger code for getting inputs of Direct Trust and other parameters from the nodes in order to calculate
the Indirect Trust.

Inputs needed : Direct Trust of the nodes, Tolerance threshold, Rate of
Increase, Rate of Decrease, Previous Similarity parameter between node A
and B (if any), Similartiy Threshold, Number of Nodes between A and B to
calculate Indirect trust of B computed by A
%}

tol_thr = 0; %% Tolerance threshold, value to be given, currently taken a 0 just for reference
r_inc = 0; %% Rate of increase, value to be given, currently taken as 0 just for reference
r_dec = 0; %% Rate of decrease, value to be given, currently taken as 0 just for reference

e_err_B_A = ((DT_x_A - DT_x_B)^2)/(DT_x_A + DT_x_B); %%Evaluation error between 2 nodes A and B, x signifes all nodes between A and B

%%Similarity parameters

s_b_a = 0; %% Similarity parameter between node A and B, value to be given, currently taken as 0 just for reference
s_b_a_prev = 0; %% previous Similarity parameter between node A and B (if any), value to be given, currently taken as 0 just for reference

if e_err_B_A < tol_thr
    s_b_a = s_b_a_prev + ((1 - s_b_a_prev)/r_inc);
else
    s_b_a = s_b_a_prev - (s_b_a_prev/r_dec);
end


%% Recommendation credibilty

sim_thr = 0; %% Similarity threshold, value to be given, currently taken as 0 just for reference
rec_b_a = 0 ; %% Recommendation credibility of Node B computed by Node A, value to be given, currently taken as 0 just for reference
if s_b_a > sim_thr
    rec_b_a = 1 - (log(s_b_a)/log(sim_thr));
else
    rec_b_a = 0;
end

%% Indirect trust (Dempster-Shaffer (D-S) theory)
%% to be calcualted based on the Node diagram, for now Let's consider N nodes between A and B
N=0; %% Value to be given, 0 is taken as ref. N is Number of nodes between A and B, if 1 node, let's name it X, it 2 then X and Y, if 3 then X,Y and Z, value to be given
if N == 1
    m_B_X = (rec_X_A * dt_B_X) / rec_X_A; %% Mass function of the Belief function of X
    it_B_A = m_B_A; %%Indirect Trust of B computed by A
elseif N==2
    m_B_X = (rec_X_A * dt_B_X) / (rec_X_A + rec_Y_A); %% Mass function of the Belief function of X
    m_B_Y = (rec_Y_A * dt_B_Y) / (rec_X_A + rec_Y_A); %% Mass function of the Belief function of Y
    it_B_A = (m_B_X + m_B_Y) - (m_B_X * m_B_Y); %%Indirect Trust of B computed by A
elseif N==3
    m_B_X = (rec_X_A * dt_B_X) / (rec_X_A + rec_Y_A + rec_Z_A); %% Mass function of the Belief function of X
    m_B_Y = (rec_Y_A * dt_B_Y) / (rec_X_A + rec_Y_A + rec_Z_A); %% Mass function of the Belief function of Y
    m_B_Z = (rec_Z_A * dt_B_Z) / (rec_X_A + rec_Y_A + rec_Z_A); %% Mass function of the Belief function of Z
    it_B_A = (m_B_X + m_B_Y + m_B_Z) - ((m_B_X * m_B_Y) + (m_B_X * m_B_Z) + (m_B_Y * m_B_Z)) + (m_B_X * m_B_Y * m_B_Z); %%Indirect Trust of B computed by A
end