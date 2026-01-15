function indx = particle_fitness_check(particles)
    % returns the index of the unfit particles
    
     indx = find(particles.Q<0);

end