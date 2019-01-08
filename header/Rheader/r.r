import click
@click.group()
def cli():
    pass

@cli.command()
def initdb():
    click.echo('Initialized the database')

@cli.command()
def dropdb():
    click.echo('Dropped the database')
#You would then invoke the Group in your setuptools entry points or other invocations:

if __name__ == '__main__':
    cli()
