# @Date:   2019-11-20T15:43:52+08:00
# @Email:  1730416009@stu.suda.edu.cn
# @Filename: CommandLine.py
# @Last modified time: 2019-11-20T17:43:44+08:00
import click
import logging


logger = logging.getLogger("RetrieveSIFTS")
streamHandler = logging.StreamHandler()
fileHandler = logging.FileHandler(filename=r"C:\Users\Nature\Downloads\test.log")
logger.setLevel(logging.DEBUG)
streamHandler.setLevel(logging.WARNING)
fileHandler.setLevel(logging.DEBUG)
formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s %(message)s")
streamHandler.setFormatter(formatter)
fileHandler.setFormatter(formatter)
logger.addHandler(streamHandler)
logger.addHandler(fileHandler)


@click.command()
@click.option("--count", default=1, help="Number of greetings.")
@click.option("--name", prompt="Your name", help="The person to greet.")
# @click.argument("name")
def hello(count, name):
    """Simple program that greets NAME for a total of COUNT times."""
    for _ in range(count):
        click.echo(click.style("Hello, %s!" % name, fg='green'))


@click.group()
def cli():
    pass


@cli.command()
@click.option("--proceed/--no-proceed", default=True, help="Whether to proceed after saving the site info.", is_flag=True)
def initdb(proceed):
    click.echo("proceed:%s" % proceed)
    click.echo('Initialized the database')


@cli.command()
def dropdb():
    click.echo('Dropped the database')
    logger.debug('This is a customer debug message')
    logger.info('This is an customer info message')
    logger.warning('This is a customer warning message')
    logger.error('This is an customer error message')
    logger.critical('This is a customer critical message')


cli.add_command(initdb)
cli.add_command(dropdb)


if __name__ == '__main__':
    # hello()
    cli()
