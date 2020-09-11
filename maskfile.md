# Tasks

Please notice this is a [mask markdown file][mask]. You need nothing special for reading it
and using the commands provided here, but if you want a bit more of automation,
then you will need to [install mask][mask-install].

If you choose to do so, then running `mask --help` will give you an overview of the tasks you can run.

## validate

> Validates the application works as expected from a functional point of view.
> There are several subcommands that you might discover by issuing `mask validate --help`.  

Source code changes are common as part of the lifecycle of whatever project. As is in the case of many
bioinformatics tools, unit and integration tests are very limited. This project is no different.

We have ran some examples provided in the [examples](examples) folder, and generated a set of checksums and
output files that you can find in the [validation](validation) folder.

We provide validation scripts to make sure that the application behaves as expected.

### all

> Run all the validations. This task might take long, but it might be the best way to validate
> the expected functionality.

~~~sh
validation/validate.sh
~~~

[mask]: https://github.com/jakedeichert/mask
[mask-install]: https://github.com/jakedeichert/mask/#installation

