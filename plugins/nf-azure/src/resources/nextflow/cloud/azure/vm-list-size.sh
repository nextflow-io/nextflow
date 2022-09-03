for region in $(az account list-locations | jq -r '.[].name'); do
  echo "Listing mv size for $region" 
  az vm list-sizes -l $region > "vm-list-size-$region.json" 
done