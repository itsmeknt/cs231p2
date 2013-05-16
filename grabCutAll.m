allImages = { %'banana1.bmp'; 'banana2.bmp'; 'banana3.bmp'; 'book.bmp'; 'bool.jpg'; 'bush.jpg'; 'ceramic.bmp'; 'cross.jpg'; 'doll.bmp'; 'elefant.bmp';
    % 'flower.jpg'; 'fullmoon.bmp'; 'grave.jpg'; 'llama.bmp'; 'memorial.jpg'; 'music.JPG'; 'person1.jpg'; 'person2.bmp'; 'person3.jpg'; 'person4.jpg';
    % 'person5.jpg'; 'person6.jpg'; 'person7.jpg'; 'person8.bmp'; 'scissors.JPG'; 'sheep.jpg'; 'stone1.JPG'; 'stone2.JPG'; 'teddy.jpg'; 'tennis.jpg'};

    'sheep.jpg', 'stone2.JPG', 'flower.jpg', 'fullmoon.bmp', 'banana1.bmp', 'cross.jpg'}; 
    

for i=1:length(allImages)
    try
        image = allImages{i}
        grabcut(image);
    catch
    end
end