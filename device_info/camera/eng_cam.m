%% Output

text_disp = {
    sprintf('Engineering Cameras');
    sprintf('\n\n     Tech Specs');
    sprintf('\n     MAIN FUNCTION:  Used for driving around on Mars and for positioning the tools on the robotic arm');
    sprintf('     LOCATION:  Various places on the rover');
    sprintf('     WEIGHT:  Less than 425 grams (less than a pound)');
    sprintf('     IMAGE SIZE:  5120 x 3840 pixels');
    sprintf('     IMAGE RESOLUTION:  20 megapixel');
    
    sprintf('\n\nHazard Avoidance Cameras (HazCams):');
    sprintf('\n\n     Tech Specs');
    sprintf('\n     MAIN FUNCTION:  Aid in autonomous navigation and obstacle avoidance');
    sprintf('     LOCATION:  Mounted at the front and rear of the rover''s body, pointing down toward the ground, about 27 inches (68 centimeters)');
    sprintf('                above ground; front: about 6.54 inches between the center of left and right eyes; back: 3.9 inches (10 centimeters),');
    sprintf('                about 31 inches (78 centimeters) above ground level');
    
    sprintf('\n\n    Mars 2020 carries six newly developed Hazard Detection Cameras, called HazCams: four on the front and two on the rear of the');
    sprintf('    rover body. HazCams detect hazards to the front and back pathways of the rover, such as large rocks, trenches, or sand dunes.');    
    sprintf('\n    Engineers also use the front HazCams to see where to move the robotic arm to take measurements, photos, and collect rock and soil samples.');
    sprintf('\n    When driving, the rover stops frequently to take new stereo images of the path ahead to evaluate potential hazards. The 3D views give Mars');
    sprintf('    2020 the ability to make its own decisions about where to drive without consulting on every move with the rover team on Earth.');
    
    
    sprintf('\n\nNavigation Cameras (NavCams):');
    sprintf('\n\n     Tech Specs');
    sprintf('\n     MAIN FUNCTION:  Aid in autonomous navigation');
    sprintf('     LOCATION:  Mounted at the front and rear of the rover''s body, pointing down toward the ground; left and right "eyes" in each');
    sprintf('                set are about 16.5 inches (42 centimeters) apart');
    
    sprintf('\n\n    Two sets of color stereo Navigation Cameras, called NavCams, help engineers navigate Mars 2020 safely, particularly when the');
    sprintf('    rover operates autonomously, making its own navigation decisions without consulting controllers on Earth.');
    
    sprintf('\n    Located up high on the rover''s mast, these two sets of black-and-white stereo cameras help engineers drive the rover around Mars.');
    sprintf('    They can see an object as small as a golf ball from 82 feet (25 meters) away. Before Mars 2020 "drives blind”, the navigation cameras initially');
    sprintf('    help ensure a safe path. Blind-drive mode occurs when engineers command the rover to drive a certain distance in a certain direction, and the');
    sprintf('    rover''s computer "brains" calculate distance from wheel rotations without looking or checking for wheel slippage.');
    
    
    sprintf('\n\nCacheCam: New Camera to Record Sample Collection');
    sprintf('\n\n     Tech Specs');
    sprintf('\n     MAIN FUNCTION:  To see down into the top of a sample tube after the sample is gathered; to take microscopic pictures of the top');
    sprintf('                     of the sample material before the tube is sealed.');
    sprintf('     LOCATION:  Inside the rover underbelly, at the top of the sample cache');
    
    sprintf('\n\n    The "CacheCam" is a single camera that looks down at the top of the sample cache. It takes pictures of sampled materials and the sample tubes');
    sprintf('    as they are being prepared for sealing and caching. This helps scientists “watch over” the samples as they are being obtained, and keeps a ');
    sprintf('    record of the entire process for each sample collected.\n\n');
};

app.TextArea.Value = text_disp;