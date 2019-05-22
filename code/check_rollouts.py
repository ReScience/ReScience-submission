from torchvision import transforms
from torchvision.utils import save_image
import torch

from utils.misc import RED_SIZE
## WARNING : THIS SHOULD BE REPLACE WITH PYTORCH 0.5
from data.loaders import RolloutObservationDataset




transform_train = transforms.Compose([
    transforms.ToPILImage(),
    transforms.Resize((RED_SIZE, RED_SIZE)),
    transforms.RandomHorizontalFlip(),
    transforms.ToTensor(),
])


dataset_train = RolloutObservationDataset('tmpexpdir',
                                          transform_train, train=True)

train_loader = torch.utils.data.DataLoader(
    dataset_train, batch_size=64, shuffle=True, num_workers=2)

for batch in train_loader:
	save_image(batch.view(64, 3, RED_SIZE, RED_SIZE), "checkrollout.png")
	break
	